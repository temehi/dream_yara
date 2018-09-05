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

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "d_misc_options.h"
#include "d_kdx_filter.h"
#include "d_bloom_filter.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString      contigsDir;
    CharString      filterFile;

    uint32_t        kmerSize;
    uint32_t        numberOfBins;
    uint64_t        bloomFilterSize;
    uint32_t        numberOfHashes;
    unsigned        threadsCount;
    bool            verbose;

    FilterType      filterType;

    std::vector<std::string> filterTypeList;

    Options() :
    kmerSize(20),
    numberOfBins(64),
    bloomFilterSize(8589934592 + filterMetadataSize), // 1GB
    numberOfHashes(4),
    threadsCount(1),
    verbose(false),
    filterType(BLOOM),
    filterTypeList({"bloom", "kmer_direct", "none"})
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
    setDefaultValue(parser, "threads", options.threadsCount);

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", options.filterTypeList);
    setDefaultValue(parser, "filter-type",  options.filterTypeList[options.filterType]);


    addOption(parser, ArgParseOption("k", "kmer-size", "The size of kmers for bloom_filter",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "kmer-size", "14");
    setMaxValue(parser, "kmer-size", "32");

    addOption(parser, ArgParseOption("nh", "num-hash", "Specify the number of hash functions to use for the bloom filter.", ArgParseOption::INTEGER));
    setMinValue(parser, "num-hash", "2");
    setMaxValue(parser, "num-hash", "5");
    setDefaultValue(parser, "num-hash", options.numberOfHashes);

    addOption(parser, ArgParseOption("bs", "bloom-size", "The size of bloom filter in GB.", ArgParseOption::INTEGER));
    setMinValue(parser, "bloom-size", "1");
    setMaxValue(parser, "bloom-size", "512");
    setDefaultValue(parser, "bloom-size", 1);
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
    getArgumentValue(options.contigsDir, parser, 0);

    // Append trailing slash if it doesn't exist.
    appendTrailingSlash(options.contigsDir);

    // Parse contigs index prefix.
    getOptionValue(options.filterFile, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.filterFile = trimExtension(options.contigsDir);
        append(options.filterFile, "bloom.filter");
    }

    getOptionValue(options.filterType, parser, "filter-type", options.filterTypeList);

    if (isSet(parser, "number-of-bins")) getOptionValue(options.numberOfBins, parser, "number-of-bins");
    if (isSet(parser, "kmer-size")) getOptionValue(options.kmerSize, parser, "kmer-size");
    if (isSet(parser, "threads")) getOptionValue(options.threadsCount, parser, "threads");
    if (isSet(parser, "num-hash")) getOptionValue(options.numberOfHashes, parser, "num-hash");

    uint64_t bloomSize;
    if (getOptionValue(bloomSize, parser, "bloom-size"))
    {
        if ((bloomSize & (bloomSize - 1)) == 0)
        {
            options.bloomFilterSize = bloomSize * 8589934592 + filterMetadataSize; // 8589934592 = 1GB
        }
        else
        {
            std::cerr <<"[ERROR] --bloom-size (-bs) parameter should be a power of 2!" << std::endl;
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
    std::string comExt = commonExtension(options.contigsDir, options.numberOfBins);

    Timer<double>       timer;
    Timer<double>       globalTimer;
    start (timer);
    start (globalTimer);

    uint32_t numThr = options.threadsCount;
    uint32_t batchSize = options.numberOfBins/numThr;
    if(batchSize * numThr < options.numberOfBins) ++batchSize;

    std::vector<std::future<void>> tasks;

    for (uint32_t taskNo = 0; taskNo < numThr; ++taskNo)
    {
        tasks.emplace_back(std::async([=, &filter] {
            for (uint32_t binNo = taskNo*batchSize; binNo < options.numberOfBins && binNo < (taskNo +1) * batchSize; ++binNo)
            {
                Timer<double>       binTimer;
                start (binTimer);

                CharString fastaFile;
                appendFileName(fastaFile, options.contigsDir, binNo);
                append(fastaFile, comExt);

                filter.addFastaFile(fastaFile, binNo);

                stop(binTimer);
                if (options.verbose)
                {
                    mtx.lock();
                    std::cerr <<"[bin " << binNo << "] Done adding kmers!\t\t\t" << binTimer << std::endl;
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
    filter.save(toCString(options.filterFile));
    stop(timer);
    if (options.verbose)
        std::cerr <<"Done saving filter (" << filter.size_mb() <<" MB)\t\t" << timer << std::endl;
    stop(globalTimer);
    std::cerr <<"\nFinshed in \t\t\t" << globalTimer << std::endl;
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
    if (!checkOutputFile(options.filterFile))
        return 1;

    try
    {
        if (options.filterType == BLOOM)
        {
            SeqAnBloomFilter<> filter  (options.numberOfBins,
                                        options.numberOfHashes,
                                        options.kmerSize,
                                        options.bloomFilterSize);

            build_filter(options, filter);
        }
        else if (options.filterType == KMER_DIRECT)
        {
            uint64_t vec_size = (1u << (2 * options.kmerSize));
            vec_size *= options.numberOfBins;
            vec_size += filterMetadataSize;
            SeqAnKDXFilter<> filter (options.numberOfBins, options.kmerSize, vec_size);

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
