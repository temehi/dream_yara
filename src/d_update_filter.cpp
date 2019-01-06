// ==========================================================================
//                                 d_update_filter.cpp
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#define UPDATE_FILTER
// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------
#include <seqan/index.h>
#include <seqan/binning_directory.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "d_misc_options.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString                                      filter_file;
    std::vector<std::pair<uint32_t, std::string> >  bin_contigs;
    uint32_t                                        number_of_bins;
    unsigned                                        threads_count;
    bool                                            verbose;
    FilterType                                      filter_type;
    std::vector<std::string>                        filter_type_list;

    Options() :
    number_of_bins(1),
    threads_count(1),
    verbose(false),
    filter_type(BLOOM),
    filter_type_list({"bloom", "kmer_direct", "none"})
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "dream_yara_update_filter");
    setShortDescription(parser, "Update Bloom Filter for DREAM-Yara");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIBLOOM-FILTER FILE \\fP> <\\fI4.fna\\fP> <\\fI7.fna\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BLOOM FILTER"));
    setValidValues(parser, 0, "filter");
    setHelpText(parser, 0, "The path of the bloom filter to be updated.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILES", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "The fasta files of the bins to updated. File names should be exactly the same us bin number (0-indexing). e.g. 0.fna");

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", options.filter_type_list);
    setDefaultValue(parser, "filter-type",  options.filter_type_list[options.filter_type]);


    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use (valid for bloom filter only).", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threads_count);

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    // Parse bloom filter path.
    getArgumentValue(options.filter_file, parser, 0);

    getOptionValue(options.filter_type, parser, "filter-type", options.filter_type_list);

    if (isSet(parser, "threads")) getOptionValue(options.threads_count, parser, "threads");

    // std::map<uint32_t, CharString>  bin_contigs;
    // Parse read input files.
    options.number_of_bins = getArgumentValueCount(parser, 1);
    options.bin_contigs.resize(options.number_of_bins);
    for (uint32_t i = 0; i < options.number_of_bins; ++i)
    {
        std::string  current_file;
        uint32_t     cur_bin_number;

        getArgumentValue(current_file, parser, 1, i);

        if (get_bin_number_from_file(cur_bin_number, current_file))
            options.bin_contigs[i] = std::make_pair(cur_bin_number, current_file);
        else
        {
            std::cerr << "File: " << current_file << "\ndoesn't have a valid name\n";
            exit(1);
        }

    }
    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function verify_fna_files()
// ----------------------------------------------------------------------------
inline bool verify_fna_files(std::vector< std::pair <uint32_t,std::string> > const & file_list)
{
    for (auto fasta_file : file_list)
    {
        if(!verify_fna_file(fasta_file.second))
            return false;
    }
    return true;
}
// ----------------------------------------------------------------------------
// Function update_filter()
// ----------------------------------------------------------------------------
template <typename TFilter>
inline void update_filter(Options & options, TFilter & filter)
{
    uint32_t no_of_bins = filter.noOfBins;

    // clear the bins to updated;
    std::vector<uint32_t> bins2update = {};
    for(int i=0; i < options.number_of_bins; ++i)
    {
        if(options.bin_contigs[i].first >= no_of_bins)
        {
            std::cerr <<"The provided bloom filter has only " << no_of_bins <<" Bins.\nRetry after removing " << options.bin_contigs[i].second << " from arguments!" << std::endl;
            exit(1);
        }
        bins2update.push_back(options.bin_contigs[i].first);
    }
    filter.clear(bins2update, options.threads_count);

    std::vector<std::future<void>> tasks;
    uint32_t batch_size = options.number_of_bins/options.threads_count;
    if(batch_size * options.threads_count < options.number_of_bins) ++batch_size;


    Timer<double>       timer;
    Timer<double>       global_timer;
    start (timer);
    start (global_timer);

    for (uint32_t task_number = 0; task_number < options.threads_count; ++task_number)
    {
        tasks.emplace_back(std::async([=, &filter, &options] {
            for (uint32_t file_number = task_number*batch_size;
                 file_number < options.number_of_bins && file_number < (task_number +1) * batch_size;
                 ++file_number)
            {
                Timer<double>       bin_timer;
                start (bin_timer);
                // read everything as CharString to avoid impure sequences crashing the program
                CharString id, seq;
                SeqFileIn seq_file_in;
                uint32_t kmer_size = filter.kmerSize;

                if (!open(seq_file_in, toCString(options.bin_contigs[file_number].second)))
                {
                    CharString msg = "Unable to open contigs file: ";
                    append(msg, CharString(options.bin_contigs[file_number].second));
                    std::cerr << msg << std::endl;
                    throw toCString(msg);
                }
                while(!atEnd(seq_file_in))
                {
                    readRecord(id, seq, seq_file_in);
                    if(length(seq) < kmer_size)
                        continue;
                    insertKmer(filter, seq, options.bin_contigs[file_number].first);
                }
                close(seq_file_in);
                stop(bin_timer);

                if (options.verbose)
                {
                    mtx.lock();
                    std::cerr   <<"[bin " << options.bin_contigs[file_number].first << "] updated using "
                                << options.bin_contigs[file_number].second << "!\t\t" << bin_timer << std::endl;
                    mtx.unlock();
                }

            }}));
    }

    for (auto &&task : tasks)
    {
        task.get();
    }

    stop(timer);
    if (options.verbose)
        std::cerr <<"All given bins are updated!\t\t" << timer << std::endl;

    start(timer);
    store(filter, toCString(options.filter_file));
    stop(timer);
    if (options.verbose)
        std::cerr <<"Done saving filter (" << size(filter) <<" MB)\t\t" << timer << std::endl;
    stop(global_timer);
    std::cerr <<"\nFinshed in \t\t\t" << global_timer << std::endl;
}
// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // verify all the new fasta files
    if (!verify_fna_files(options.bin_contigs))
        return 1;

    try
    {
        if (options.filter_type == BLOOM)
        {
            typedef BinningDirectory<InterleavedBloomFilter,
                                     BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;
            BinningDirectoriesIBF filter;
            retrieve(filter, toCString(options.filter_file));
            update_filter(options, filter);
        }
        else if (options.filter_type == KMER_DIRECT)
        {
            typedef BinningDirectory<DirectAddressing,
                                     BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesDA;

            BinningDirectoriesDA filter;
            retrieve(filter, toCString(options.filter_file));
            update_filter(options, filter);
        }
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
