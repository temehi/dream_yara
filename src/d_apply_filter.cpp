// ==========================================================================
//      DREAM-Yara - Yet Another Read Aligner within DREAM framework
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

#define YARA_APPLY_FILTER

// ============================================================================
// Forwards
// ============================================================================

struct Options;

// ============================================================================
// Prerequisites
// ============================================================================
//#include "BloomFilter.hpp"
//#include "ntHashIterator.hpp"

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <vector>
#include <string>
#include <random>
#include <math.h>
#include <fstream>
#include <iostream>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>
#include <seqan/binning_directory.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "basic_alphabet.h"
#include "file_pair.h"
#include "file_prefetched.h"
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "index_fm.h"
#include "bits_reads.h"
#include "bits_hits.h"
#include "bits_context.h"
#include "bits_matches.h"
#include "bits_seeds.h"
#include "bits_bucket.h"
#include "find_verifier.h"
#include "find_extender.h"
#include "misc_options.h"
#include "d_misc_options.h"
#include "mapper_collector.h"
#include "mapper_classifier.h"
#include "mapper_ranker.h"
#include "mapper_filter.h"
#include "mapper_extender.h"
#include "mapper_verifier.h"
#include "mapper_aligner.h"
#include "mapper_writer.h"
#include "mapper.h"

#include "d_apply_filter.h"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, FilterAppOptions const & options)
{
    setAppName(parser, "filter_app");
    setShortDescription(parser, "Apply Filter to a set of reads.");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    // Setup mandatory arguments.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX DIRECTORY\\fP> <\\fISE-READS FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FILTER FILE"));
    setHelpText(parser, 0, "The path to a bloom filter.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS FILE", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one single-end or two paired-end / mate-pair read files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays global statistics."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Displays extensive statistics for each batch of reads."));

    // Setup output options.
    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.",
                                     ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", "csv");

    addOption(parser, ArgParseOption("e", "error-rate", "Consider alignments within this percentual number of errors. \
                                     Increase this threshold to increase the number of mapped reads. \
                                     Decrease this threshold to decrease the runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", 100.0 * options.error_rate);

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threads_count);

    addOption(parser, ArgParseOption("rb", "reads-batch", "Specify the number of reads to process in one batch.",
                                     ArgParseOption::INTEGER));

    setMinValue(parser, "reads-batch", "1000");
    setMaxValue(parser, "reads-batch", "5000000");
    setDefaultValue(parser, "reads-batch", options.reads_count);
    hideOption(getOption(parser, "reads-batch"));

    addOption(parser, ArgParseOption("s", "stat-only", "Do not write output to file. Print out stat only! output-file will be ignored!"));

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", options.filter_typeList);
    setDefaultValue(parser, "filter-type",  options.filter_typeList[options.filter_type]);
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(FilterAppOptions & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse the filter file
    getArgumentValue(options.filter_file, parser, 0);


    getArgumentValue(options.reads_file.i1, parser, 1, 0);

    // Parse output file.
    getOptionValue(options.output_file, parser, "output-file");

    // Parse mapping options.
    unsigned error_rate;
    if (getOptionValue(error_rate, parser, "error-rate"))
        options.error_rate = error_rate / 100.0;

    // Parse performance options.
    getOptionValue(options.threads_count, parser, "threads");
    getOptionValue(options.reads_count, parser, "reads-batch");

    getOptionValue(options.filter_type, parser, "filter-type", options.filter_typeList);
    getOptionValue(options.stat_only, parser, "stat-only");

    if (isSet(parser, "verbose")) options.verbose = 1;
    if (isSet(parser, "very-verbose")) options.verbose = 2;

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function check_read_files()
// ----------------------------------------------------------------------------
//
bool check_read_files(FilterAppOptions const &  options)
{
    // check if read file(s) exist(s)
    SeqFileIn seqFile;
    if (!open(seqFile, toCString(options.reads_file.i1)))
    {
        std::cerr << "Unable to open read file "<< toCString(options.reads_file.i1) <<"!\n";
        return false;
    }
    close(seqFile);
    return true;
}

// ----------------------------------------------------------------------------
// Function read_filter_metadata()
// ----------------------------------------------------------------------------
bool read_filter_metadata(FilterAppOptions &  options)
{
    std::ifstream in(toCString(options.filter_file), std::ios::in | std::ios::binary);
    uint64_t x = BD_METADATA_SIZE/8; //bits -> bytes

    sdsl::int_vector<64>  metadata_vec(BD_METADATA_SIZE/64, 0); //bits -> uint64_t
    in.seekg(-x, in.end); // seek from end of file

    uint64_t* p  = &(metadata_vec[0]);
    in.read((char*)p, x * sizeof(uint64_t));

    //    std::cout << metadata_vec << std::endl;
    options.number_of_bins = metadata_vec[0];
    options.kmer_size = metadata_vec[2];

    return true;
}

// ----------------------------------------------------------------------------
// Function init_filter_app()
// ----------------------------------------------------------------------------
template <typename TThreading>
inline void init_filter_app(FilterAppOptions & options, TThreading const & threading)
{
    typedef ReadMapperConfig<TThreading>  TConfig;
    FilterApp<void, TConfig> me(options);

    start(me.timer);

    if (me.options.filter_type == BLOOM)
    {
        typedef BinningDirectory<InterleavedBloomFilter,
        BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;
        BinningDirectoriesIBF filter;
        retrieve(filter, toCString(me.options.filter_file));

        me.options.kmer_size = filter.kmerSize;
        me.options.number_of_bins = filter.noOfBins;

        stop(me.timer);
        me.options.load_filter_time += getValue(me.timer);
        run_filter_app(me, filter);
    }
    else if (me.options.filter_type == KMER_DIRECT)
    {
        typedef BinningDirectory<DirectAddressing,
        BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesDA;
        BinningDirectoriesDA filter;
        retrieve(filter, toCString(me.options.filter_file));

        me.options.kmer_size = filter.kmerSize;
        me.options.number_of_bins = filter.noOfBins;

        stop(me.timer);
        me.options.load_filter_time += getValue(me.timer);
        run_filter_app(me, filter);
    }
    else
    {
        // dummy filter in case of nofilter option
        typedef BinningDirectory<InterleavedBloomFilter,
        BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;

        BinningDirectoriesIBF filter(64, 3, 20, 1);

        stop(me.timer);
        me.options.load_filter_time += getValue(me.timer);
        run_filter_app(me, filter);
    }
    std::cerr << "Number of reads:\t\t" << (double)me.options.reads_count << std::endl;
    std::cerr << "Avg reads per bin:\t\t" << (double)me.options.passed_reads_count / me.options.number_of_bins << std::endl;
    std::cerr << "Filter loading time:\t" << me.options.load_filter_time << " sec" << std::endl;
    std::cerr << "Filter reads time:\t\t" << me.options.filter_reads_time << " sec" << std::endl;
}


// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
//
int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    FilterAppOptions options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;


    if (!check_read_files(options))
        return 1;

    try
    {

#ifdef _OPENMP
    if (options.threadsCount > 1)
        init_filter_app(options, Parallel());
    else
#endif
        init_filter_app(options, Serial());
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
