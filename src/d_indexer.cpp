// ==========================================================================
//                                 d_indexer.cpp
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


#define YARA_INDEXER
// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "d_misc_options.h"
#include "index_fm.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------
struct Options
{
    std::vector<std::pair<uint32_t, std::string> >  bin_contigs;
    CharString                                      contigsIndexFile;
    CharString                                      cur_index_file;

    uint32_t                                        number_of_bins;
    uint32_t                                        cur_bin_number;

    uint64_t                                        contigsSize;
    uint64_t                                        contigsMaxLength;
    uint64_t                                        contigsSum;

    unsigned                                        threads_count;

    bool                                            verbose;

    Options() :
    number_of_bins(1),
    cur_bin_number(0),
    contigsSize(),
    contigsMaxLength(),
    contigsSum(),
    threads_count(1),
    verbose(false)
    {}
};

// ----------------------------------------------------------------------------
// Class YaraIndexer
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = void>
struct YaraIndexer
{
    typedef SeqStore<TSpec, YaraContigsConfig<> >   TContigs;

    Options const &     options;
    TContigs            contigs;
    SeqFileIn           contigsDir;
    Timer<double>       timer;

    YaraIndexer(Options const & options) :
    options(options)
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
    setAppName(parser, "dream_yara_indexer");
    setShortDescription(parser, "DREAM-Yara Indexer");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fI0.fna\\fP> <\\fI1.fna\\fP>");

//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));
//    //    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
//    setHelpText(parser, 0, "A directory containing reference genome files.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILES", true));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "The fasta files of the bins to updated. File names should be exactly the same us bin number (0-indexing). e.g. 0.fna");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-prefix", "Specify a filename prefix for the reference genome index. \
                                     Default: use the filename prefix of the reference genome.", ArgParseOption::OUTPUT_PREFIX));
    setRequired(parser, "output-prefix");

    addOption(parser, ArgParseOption("td", "tmp-dir", "Specify a temporary directory where to construct the index. \
                                     Default: use the output directory.", ArgParseOption::STRING));
//    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
//                                     ArgParseOption::INTEGER));
//    setMinValue(parser, "number-of-bins", "1");
//    setMaxValue(parser, "number-of-bins", "10000");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use (valid for bloom filter only).", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "32");
    setDefaultValue(parser, "threads", options.threads_count);
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
//    getArgumentValue(options.contigsDir, parser, 0);

    // std::map<uint32_t, CharString>  binContigs;
    // Parse contig input files.
    options.number_of_bins = getArgumentValueCount(parser, 0);
    options.bin_contigs.resize(options.number_of_bins);
    for (uint32_t i = 0; i < options.number_of_bins; ++i)
    {
        std::string  current_file;
        uint32_t     cur_bin_number;

        getArgumentValue(current_file, parser, 0, i);

        if (get_bin_number_from_file(cur_bin_number, current_file))
            options.bin_contigs[i] = std::make_pair(cur_bin_number, current_file);
        else
        {
            std::cerr << "File: " << current_file << "\ndoesn't have a valid name\n";
            exit(1);
        }

    }

    // Parse contigs index prefix.
    getOptionValue(options.contigsIndexFile, parser, "output-prefix");

    // Parse and set temp dir.
    CharString tmp_dir;
    getOptionValue(tmp_dir, parser, "tmp-dir");
    if (!isSet(parser, "tmp-dir"))
    {
        tmp_dir = getPath(options.contigsIndexFile);
        if (empty(tmp_dir))
            getCwd(tmp_dir);
    }
    setEnv("TMPDIR", tmp_dir);

//    if (isSet(parser, "number-of-bins")) getOptionValue(options.number_of_bins, parser, "number-of-bins");
    if (isSet(parser, "threads")) getOptionValue(options.threads_count, parser, "threads");

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function load_contigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void load_contigs(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.verbose)
    {
        mtx.lock();
        std::cerr << "[bin " << me.options.cur_bin_number << "] Loading reference ..." << std::endl;
        mtx.unlock();
    }

    if (!open(me.contigsDir, toCString(me.options.cur_index_file)))
        throw RuntimeError("Error while opening the reference file.");

    try
    {
        readRecords(me.contigs, me.contigsDir);
        trimSeqNames(me.contigs);
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }

}

// ----------------------------------------------------------------------------
// Function save_contigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void save_contigs(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.verbose)
    {
        mtx.lock();
        std::cerr <<"[bin " << me.options.cur_bin_number << "] Saving reference ..." << std::endl;
        mtx.unlock();
    }

    if (!saveContigsLimits(me.options) || !save(me.contigs, toCString(me.options.contigsIndexFile)))
        throw RuntimeError("Error while saving the reference.");
}

// ----------------------------------------------------------------------------
// Function save_index()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TContigsSum, typename TSpec, typename TConfig>
void save_index(YaraIndexer<TSpec, TConfig> & me)
{
    typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum>    TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;

    if (me.options.verbose)
    {
        mtx.lock();
        std::cerr << "[bin " << me.options.cur_bin_number << "] Building reference index ..." << std::endl;
        mtx.unlock();
    }


    // Randomly replace Ns with A, C, G, T.
    randomizeNs(me.contigs);

    // IndexFM is built on the reversed contigs.
    reverse(me.contigs);

    TIndex index;

    // This assignment *copies* the contigs to the index as the types differ.
    setValue(index.text, me.contigs.seqs);

    // Clear the contigs - the index now owns its own copy.
    clear(me.contigs);
    shrinkToFit(me.contigs);

    try
    {
        // Iterator instantiation triggers index construction.
        typename Iterator<TIndex, TopDown<> >::Type it(index);
        ignoreUnusedVariableWarning(it);
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to index the reference.");
    }
    catch (IOError const & /* e */)
    //    catch (PageFrameError const & /* e */)
    {
        throw RuntimeError("Insufficient disk space to index the reference. \
                           Specify a bigger temporary folder using the options --tmp-dir.");
    }

    if (me.options.verbose)
    {
        mtx.lock();
        std::cerr << "[bin " << me.options.cur_bin_number << "] Saving reference index:\t\t\t" << std::endl;
        mtx.unlock();
    }
    if (!save(index, toCString(me.options.contigsIndexFile)))
        throw RuntimeError("Error while saving the reference index file.");
}

template <typename TContigsSize, typename TContigsLen, typename TSpec, typename TConfig>
void save_index(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        save_index<TContigsSize, TContigsLen, uint32_t>(me);
    }
    else
    {
        save_index<TContigsSize, TContigsLen, uint64_t>(me);
    }
}

template <typename TContigsSize, typename TSpec, typename TConfig>
void save_index(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        save_index<TContigsSize, uint32_t>(me);
    }
    else
    {
#ifdef DDR_YARA_LARGE_CONTIGS
        save_index<TContigsSize, uint64_t>(me);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TSpec, typename TConfig>
void save_index(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        save_index<uint8_t>(me);
    }
    else if (me.options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        save_index<uint16_t>(me);
    }
    else
    {
#ifdef DDR_YARA_LARGE_CONTIGS
        save_index<uint32_t>(me);
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

// ----------------------------------------------------------------------------
// Function runYaraIndexer()
// ----------------------------------------------------------------------------
void runYaraIndexer(Options & options)
{
    YaraIndexer<> indexer(options);
    load_contigs(indexer);
    setContigsLimits(options, indexer.contigs.seqs);
    save_contigs(indexer);
    save_index(indexer);
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

    try
    {
        std::vector<std::future<void>> tasks;
        uint32_t batch_size = options.number_of_bins/options.threads_count;
        if(batch_size * options.threads_count < options.number_of_bins) ++batch_size;

        Timer<double>       timer;
        Timer<double>       global_timer;
        start (timer);
        start (global_timer);

        for (uint32_t task_number = 0; task_number < options.threads_count; ++task_number)
        {
            tasks.emplace_back(std::async([=] {
                for (uint32_t file_number = task_number*batch_size;
                     file_number < options.number_of_bins && file_number < (task_number +1) * batch_size;
                     ++file_number)
                {
                    Timer<double>       bin_timer;
                    start (bin_timer);

                    Options bin_options = options;
                    append_file_name(bin_options.contigsIndexFile, options.contigsIndexFile, options.bin_contigs[file_number].first);
                    bin_options.cur_bin_number = options.bin_contigs[file_number].first;
                    bin_options.cur_index_file = options.bin_contigs[file_number].second;
                    runYaraIndexer(bin_options);

                    stop(bin_timer);

                    if (options.verbose)
                    {
                        mtx.lock();
                        std::cerr <<"[bin " << options.bin_contigs[file_number].first << "] Done indexing reference\t\t\t" << bin_timer << std::endl;
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
            std::cerr <<"All bins are done indexing reference!\t" << timer << std::endl;

        stop(global_timer);
        std::cerr <<"\nFinshed in \t\t\t" << global_timer << std::endl;
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
