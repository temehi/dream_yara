// ==========================================================================
//                                 d_build_filter.cpp
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

#define BUILD_FILTER
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
    CharString      contigs_dir;
    CharString      filter_file;

    uint32_t        kmer_size;
    uint32_t        number_of_bins;
    uint64_t        size_of_ibf;
    uint32_t        number_of_hashes;
    unsigned        threads_count;
    bool            verbose;

    FilterType      filter_type;

    std::vector<std::string> filter_type_list;

    Options() :
    kmer_size(20),
    number_of_bins(64),
    size_of_ibf(1_g), // 1000_m
    number_of_hashes(4),
    threads_count(1),
    verbose(false),
    filter_type(BLOOM),
    filter_type_list({"bloom", "kmer_direct", "none"})
    {
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "dream_yara_build_filter");
    setShortDescription(parser, "Build Filter for DREAM-Yara");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILES DIR \\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));
    //    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A directory containing reference genome files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output filename for the filter. \
                                     Default: use the directory name of reference genomes.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", "filter");

    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
                                     ArgParseOption::INTEGER));

    setMinValue(parser, "number-of-bins", "1");
    setMaxValue(parser, "number-of-bins", "4194300");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threads_count);

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", options.filter_type_list);
    setDefaultValue(parser, "filter-type",  options.filter_type_list[options.filter_type]);


    addOption(parser, ArgParseOption("k", "kmer-size", "The size of kmers for bloom_filter",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "kmer-size", "14");
    setMaxValue(parser, "kmer-size", "32");

    addOption(parser, ArgParseOption("nh", "num-hash", "Specify the number of hash functions to use for the bloom filter.", ArgParseOption::INTEGER));
    setMinValue(parser, "num-hash", "2");
    setMaxValue(parser, "num-hash", "5");
    setDefaultValue(parser, "num-hash", options.number_of_hashes);

    addOption(parser, ArgParseOption("bs", "bloom-size",
            "The size of bloom filter suffixed by either M or G for mega bytes or giga bytes respectively.",
            ArgParseOption::STRING));
    setDefaultValue(parser, "bloom-size", "1G");
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

    // Parse contigs input file.
    getArgumentValue(options.contigs_dir, parser, 0);

    // Append trailing slash if it doesn't exist.
    append_trailing_slash(options.contigs_dir);

    // Parse contigs index prefix.
    getOptionValue(options.filter_file, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.filter_file = trimExtension(options.contigs_dir);
        append(options.filter_file, "bloom.filter");
    }

    getOptionValue(options.filter_type, parser, "filter-type", options.filter_type_list);

    if (isSet(parser, "number-of-bins")) getOptionValue(options.number_of_bins, parser, "number-of-bins");
    if (isSet(parser, "kmer-size")) getOptionValue(options.kmer_size, parser, "kmer-size");
    if (isSet(parser, "threads")) getOptionValue(options.threads_count, parser, "threads");
    if (isSet(parser, "num-hash")) getOptionValue(options.number_of_hashes, parser, "num-hash");

    std::string ibf_size;
    if (getOptionValue(ibf_size, parser, "bloom-size"))
    {
        uint64_t base = std::stoi(ibf_size);
        switch (ibf_size.at(ibf_size.size()-1))
        {
            case 'G': case 'g':
                options.size_of_ibf =  base * 8*1024*1024*1024;
                break;
            case 'M': case 'm':
                options.size_of_ibf =  base * 8*1024*1024;
                break;
            default:
                std::cerr <<"[ERROR] invalid --bloom-size (-bs) parameter provided. (eg 256M, 1g)" << std::endl;
                exit(1);
        }

    }
    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function build_filter()
// ----------------------------------------------------------------------------
template <typename TFilter>
inline void build_filter(Options & options, TFilter & filter)
{
    std::string com_ext = common_ext(options.contigs_dir, options.number_of_bins);

    Timer<double>       timer;
    Timer<double>       global_timer;
    start (timer);
    start (global_timer);

    uint32_t batch_size = options.number_of_bins/options.threads_count;
    if(batch_size * options.threads_count < options.number_of_bins) ++batch_size;

    std::vector<std::future<void>> tasks;

    for (uint32_t task_number = 0; task_number < options.threads_count; ++task_number)
    {
        tasks.emplace_back(std::async([=, &filter] {
            for (uint32_t bin_number = task_number*batch_size;
                bin_number < options.number_of_bins && bin_number < (task_number +1) * batch_size;
                ++bin_number)
            {
                Timer<double>       bin_timer;
                start (bin_timer);

                CharString seq_file_path;
                append_file_name(seq_file_path, options.contigs_dir, bin_number);
                append(seq_file_path, com_ext);

                // read everything as CharString to avoid impure sequences crashing the program
                CharString id, seq;
                SeqFileIn seq_file_in;
                if (!open(seq_file_in, toCString(seq_file_path)))
                {
                    CharString msg = "Unable to open contigs file: ";
                    append(msg, CharString(seq_file_path));
                    std::cerr << msg << std::endl;
                    throw toCString(msg);
                }
                while(!atEnd(seq_file_in))
                {
                    readRecord(id, seq, seq_file_in);
                    if(length(seq) < options.kmer_size)
                        continue;
                    insertKmer(filter, seq, bin_number);
                }
                close(seq_file_in);
                stop(bin_timer);
                if (options.verbose)
                {
                    mtx.lock();
                    std::cerr <<"[bin " << bin_number << "] Done adding kmers!\t\t\t" << bin_timer << std::endl;
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
        std::cerr <<"All bins are done adding kmers!\t\t" << timer << std::endl;

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

    // check if file already exists or can be created
    if (!check_output_file(options.filter_file))
        return 1;


    try
    {
        if (options.filter_type == BLOOM)
        {

            typedef BinningDirectory<InterleavedBloomFilter,
                                    BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;
            BinningDirectoriesIBF filter(options.number_of_bins,
                                         options.number_of_hashes,
                                         options.kmer_size,
                                         options.size_of_ibf);
            build_filter(options, filter);
        }
         else if (options.filter_type == KMER_DIRECT)
         {
            typedef BinningDirectory<DirectAddressing,
                                    BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesDA;
             BinningDirectoriesDA filter( options.number_of_bins, options.kmer_size);
             build_filter(options, filter);
         }

    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
