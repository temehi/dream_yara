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

#define YARA_MAPPER

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

#include "d_mapper.h"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, DisOptions const & d_options)
{
    setAppName(parser, "dream_yara_mapper");
    setShortDescription(parser, "DREAM-Yara Mapper");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    // Setup mandatory arguments.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX DIRECTORY\\fP> <\\fISE-READS FILE\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX DIRECTORY\\fP> <\\fIPE-READS FILE 1\\fP> <\\fIPE-READS FILE 2\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE INDEX DIRECTORY"));
    setHelpText(parser, 0, "A directory containing multiple indices of reference genomes.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS FILE", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one single-end or two paired-end / mate-pair read files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays global statistics."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Displays extensive statistics for each batch of reads."));

    // Setup output d_options.
    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.",
                                     ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", BamFileOut::getFileExtensions());

    addOption(parser, ArgParseOption("f", "output-format", "Specify an output format. Note: when specifying the option \
                                     --output-file, the output format is taken from the filename \
                                     extension.", ArgParseOption::STRING));
    setValidValues(parser, "output-format", getExtensionsWithoutLeadingDot(BamFileOut::getFileExtensions()));
    setDefaultValue(parser, "output-format", "sam");

#if SEQAN_HAS_ZLIB
    addOption(parser, ArgParseOption("u", "uncompressed-bam", "Turn off compression of BAM written to standard output."));
    hideOption(getOption(parser, "uncompressed-bam"));
#endif

    addOption(parser, ArgParseOption("rg", "read-group", "Specify a read group for all records in the SAM/BAM file.",
                                     ArgParseOption::STRING));
    setDefaultValue(parser, "read-group", d_options.readGroup);
    addOption(parser, ArgParseOption("sm", "secondary-matches", "Specify whether to output secondary matches in \
                                     the XA tag of the primary alignment, as separate \
                                     secondary records, or to omit them.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "secondary-matches", d_options.secondaryMatchesList);
    setDefaultValue(parser, "secondary-matches", d_options.secondaryMatchesList[d_options.secondaryMatches]);

    // Keep legacy option to display an error
    addOption(parser, ArgParseOption("sa", "secondary-alignments", "This option has been renamed to 'secondary-matches'.",
                                     ArgParseOption::STRING));
    hideOption(parser, "sa");

    // Turn off for DREAM-Yara temporarly till matches and cigar managment is handeled properly
    // addOption(parser, ArgParseOption("as", "align-secondary", "Compute and output co- and suboptimal \
    //                                  match alignments. Ignored if '-sa omit' is set."));

    addOption(parser, ArgParseOption("ra", "rabema-alignments", "Output alignments compatible with the \
                                     Read Alignment BEnchMArk (RABEMA)."));

    addOption(parser, ArgParseOption("sk", "skip-sam-headers", "Skip writing SQ: headers to SAM output (works only with SAM format)."));

    // Setup mapping d_options.
    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("e", "error-rate", "Consider alignments within this percentual number of errors. \
                                     Increase this threshold to increase the number of mapped reads. \
                                     Decrease this threshold to decrease the runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", 100.0 * d_options.errorRate);

    addOption(parser, ArgParseOption("s", "strata-rate", "Consider suboptimal alignments within this percentual number \
                                     of errors from the optimal alignment. Increase this threshold to increase \
                                     the number of alternative alignments at the expense of runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "strata-rate", "0");
    setMaxValue(parser, "strata-rate", "10");
    setDefaultValue(parser, "strata-rate", 100.0 * d_options.strataRate);

    addOption(parser, ArgParseOption("y", "sensitivity", "Sensitivity with respect to edit distance. \
                                     Full sensitivity guarantees to find all considered alignments \
                                     but increases runtime, low sensitivity decreases runtime by \
                                     breaking such guarantee.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "sensitivity", d_options.sensitivityList);
    setDefaultValue(parser, "sensitivity", d_options.sensitivityList[d_options.sensitivity]);

    // Setup paired-end mapping d_options.
    addSection(parser, "Paired-End Mapping Options");

    addOption(parser, ArgParseOption("ll", "library-length", "Expected library length. Default: autodetected.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");

    addOption(parser, ArgParseOption("ld", "library-deviation", "Deviation from the expected library length. \
                                     Default: autodetected.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-deviation", "1");

    addOption(parser, ArgParseOption("i", "indel-rate", "Rescue unaligned ends within this percentual number of indels.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "indel-rate", "0");
    setMaxValue(parser, "indel-rate", "50");
    setDefaultValue(parser, "indel-rate", 100.0 * d_options.indelRate);

    addOption(parser, ArgParseOption("ni", "no-indels", "Turn off the rescue of unaligned ends containing indels."));

    //    addOption(parser, ArgParseOption("lo", "library-orientation", "Expected orientation of the segments in the library.",
    //                                     ArgParseOption::STRING));
    //    setValidValues(parser, "library-orientation", d_options.libraryOrientationList);
    //    setDefaultValue(parser, "library-orientation", d_options.libraryOrientationList[d_options.libraryOrientation]);

    //    addOption(parser, ArgParseOption("la", "anchor", "Anchor one read and verify its mate."));

    // Setup performance d_options.
    addSection(parser, "Performance Options");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", d_options.threadsCount);

    addOption(parser, ArgParseOption("rb", "reads-batch", "Specify the number of reads to process in one batch.",
                                     ArgParseOption::INTEGER));

    setMinValue(parser, "reads-batch", "1000");
    setMaxValue(parser, "reads-batch", "5000000");
    setDefaultValue(parser, "reads-batch", d_options.readsCount);
    hideOption(getOption(parser, "reads-batch"));

    // Setup Distributed mapper d_options.
    addSection(parser, "Distributed mapper Options");
//    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
//                                     ArgParseOption::INTEGER));
//    setMinValue(parser, "number-of-bins", "1");
//    setMaxValue(parser, "number-of-bins", "1024");
//    setDefaultValue(parser, "number-of-bins", d_options.numberOfBins);

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", d_options.filter_typeList);
    setDefaultValue(parser, "filter-type",  d_options.filter_typeList[d_options.filter_type]);


    addOption(parser, ArgParseOption("fi", "bloom-filter", "The path to a bloom filter. Default: will look for bloom.filter file inside the indices directory.", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "bloom-filter", "filter");
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(DisOptions & d_options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse indexed genome input file.
    getArgumentValue(d_options.indices_dir, parser, 0);

    // Append trailing slash if it doesn't exist.
    append_trailing_slash(d_options.indices_dir);

    // Parse read input files.
    switch (getArgumentValueCount(parser, 1))
    {
        case 1:
            getArgumentValue(d_options.readsFile.i1, parser, 1, 0);
            d_options.singleEnd = true;
            break;
        case 2:
            getArgumentValue(d_options.readsFile.i1, parser, 1, 0);
            getArgumentValue(d_options.readsFile.i2, parser, 1, 1);
            d_options.singleEnd = false;
            break;
        default:
            std::cerr << getAppName(parser) << ": Too many arguments!" << std::endl;
            return ArgumentParser::PARSE_ERROR;
    }

    // Parse output file.
    getOptionValue(d_options.super_output_file, parser, "output-file");

    // Parse output format.
    CharString outputFormat;
    if (getOptionValue(outputFormat, parser, "output-format"))
    {
        addLeadingDot(outputFormat);
        guessFormatFromFilename(outputFormat, d_options.outputFormat);
    }
    else
        assign(d_options.outputFormat, Sam());

#if SEQAN_HAS_ZLIB
    getOptionValue(d_options.uncompressedBam, parser, "uncompressed-bam");
#endif

    // Parse output d_options.
    getOptionValue(d_options.readGroup, parser, "read-group");
    getOptionValue(d_options.secondaryMatches, parser, "secondary-alignments", d_options.secondaryMatchesList);
    getOptionValue(d_options.rabema, parser, "rabema-alignments");

    if (isSet(parser, "secondary-alignments"))
    {
        std::cerr << getAppName(parser) << ": The 'secondary-alignments' option has been renamed to 'secondary-matches'!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Turn off for DREAM-Yara temporarly till matches and cigar managment is handeled properly
    // getOptionValue(d_options.alignSecondary, parser, "align-secondary");
    // if (d_options.alignSecondary && d_options.secondaryMatches == OMIT)
    // {
    //     d_options.alignSecondary = false;
    //     std::cerr << getAppName(parser) << ": WARNING, ignoring '-as' as '-sa omit' is set." << std::endl;
    // }

    if (isSet(parser, "skip-sam-headers")) d_options.skip_sam_header = true;

    // Parse mapping d_options.
    unsigned errorRate;
    if (getOptionValue(errorRate, parser, "error-rate"))
        d_options.errorRate = errorRate / 100.0;

    unsigned strataRate;
    if (getOptionValue(strataRate, parser, "strata-rate"))
        d_options.strataRate = strataRate / 100.0;

    getOptionValue(d_options.sensitivity, parser, "sensitivity", d_options.sensitivityList);

    // Parse paired-end mapping d_options.
    getOptionValue(d_options.libraryLength, parser, "library-length");
    getOptionValue(d_options.libraryDev, parser, "library-deviation");
    //    getOptionValue(d_options.libraryOrientation, parser, "library-orientation", d_options.libraryOrientationList);

    unsigned indelRate;
    if (getOptionValue(indelRate, parser, "indel-rate"))
        d_options.indelRate = indelRate / 100.0;

    d_options.verifyMatches = !isSet(parser, "no-indels");

    // Parse performance d_options.
    getOptionValue(d_options.threadsCount, parser, "threads");
    getOptionValue(d_options.readsCount, parser, "reads-batch");

//    // Parse Distributed mapper options
//    getOptionValue(d_options.numberOfBins, parser, "number-of-bins");
//    if (isSet(parser, "number-of-bins")) getOptionValue(d_options.numberOfBins, parser, "number-of-bins");

    // Parse contigs index prefix.
    getOptionValue(d_options.filter_file, parser, "bloom-filter");
    if (!isSet(parser, "bloom-filter"))
    {
        d_options.filter_file = d_options.indices_dir;
        append(d_options.filter_file, "bloom.filter");
    }

    getOptionValue(d_options.filter_type, parser, "filter-type", d_options.filter_typeList);

    if (isSet(parser, "verbose")) d_options.verbose = 1;
    if (isSet(parser, "very-verbose")) d_options.verbose = 2;

    // Get version.
    d_options.version = getVersion(parser);

    // Get command line.
    for (int i = 0; i < argc; i++)
    {
        append(d_options.commandLine, argv[i]);
        appendValue(d_options.commandLine, ' ');
    }
    eraseBack(d_options.commandLine);

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configure_d_mapper(DisOptions & d_options,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    if (d_options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawn_d_mapper<TContigsSize, TContigsLen, uint32_t>(d_options, threading, sequencing, distance);
    }
    else
    {
        spawn_d_mapper<TContigsSize, TContigsLen, uint64_t>(d_options, threading, sequencing, distance);
    }
}

template <typename TContigsSize, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configure_d_mapper(DisOptions & d_options,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    if (d_options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configure_d_mapper<TContigsSize, uint32_t>(d_options, threading, sequencing, distance);
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configure_d_mapper<TContigsSize, uint64_t>(d_options, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TThreading, typename TSequencing, typename TSeedsDistance>
void configure_d_mapper(DisOptions & d_options,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    d_options.contigsMaxLength = 0;
    d_options.contigsSize = 0;
    d_options.contigsSum = 0;
    d_options.contig_offsets.resize(d_options.number_of_bins, 0);
    // We aggregate individual limit here to configure the dis_mapper limits
    for (uint32_t i=0; i < d_options.number_of_bins; ++i)
    {
        d_options.contig_offsets[i] = d_options.contigsSize;
        Options options = d_options;
        append_file_name(options.contigsIndexFile, d_options.indices_dir, i);
        if (!openContigsLimits(options))
            throw RuntimeError("Error while opening contig limits file.");

       d_options.contigsMaxLength   = std::max(options.contigsMaxLength, d_options.contigsMaxLength);
       d_options.contigsSize       += options.contigsSize;
       d_options.contigsSum        += options.contigsSum;
    }

    if (d_options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configure_d_mapper<uint8_t>(d_options, threading, sequencing, distance);
    }
    else if (d_options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configure_d_mapper<uint16_t>(d_options, threading, sequencing, distance);
    }
    else
    {
        configure_d_mapper<uint32_t>(d_options, threading, sequencing, distance);
    }
}

template <typename TThreading, typename TSequencing>
void configure_d_mapper(DisOptions & d_options,
                        TThreading const & threading,
                        TSequencing const & sequencing)
{
    if (d_options.sensitivity == FULL)
        return configure_d_mapper(d_options, threading, sequencing, EditDistance());
    else
        return configure_d_mapper(d_options, threading, sequencing, HammingDistance());
}

template <typename TThreading>
void configure_d_mapper(DisOptions & d_options,
                        TThreading const & threading)
{
    if (d_options.singleEnd)
        configure_d_mapper(d_options, threading, SingleEnd());
    else
        configure_d_mapper(d_options, threading, PairedEnd());
}

void configure_d_mapper(DisOptions & d_options)
{
#ifdef _OPENMP
    if (d_options.threadsCount > 1)
        configure_d_mapper(d_options, Parallel());
    else
#endif
        configure_d_mapper(d_options, Serial());
}

// ----------------------------------------------------------------------------
// Function check_read_files()
// ----------------------------------------------------------------------------
//
bool check_read_files(DisOptions const &  d_options)
{
    // check if read file(s) exist(s)
    SeqFileIn seq_file_in;
    if (!open(seq_file_in, toCString(d_options.readsFile.i1)))
    {
        std::cerr << "Unable to open read file "<< toCString(d_options.readsFile.i1) <<"!\n";
        return false;
    }
    close(seq_file_in);
    if (!d_options.singleEnd && !open(seq_file_in, toCString(d_options.readsFile.i2)))
    {
        std::cerr << "Unable to open read file "<< toCString(d_options.readsFile.i2) <<"!\n";
        return false;
    }
    close(seq_file_in);
    return true;
}

// ----------------------------------------------------------------------------
// Function read_filter_metadata()
// ----------------------------------------------------------------------------

bool read_filter_metadata(DisOptions &  d_options)
{
    std::ifstream in(toCString(d_options.filter_file), std::ios::in | std::ios::binary);
    uint64_t x = BD_METADATA_SIZE/8; //bits -> bytes

    sdsl::int_vector<64>  metadata_vec(BD_METADATA_SIZE/64, 0); //bits -> uint64_t
    in.seekg(-x, in.end); // seek from end of file

    uint64_t* p  = &(metadata_vec[0]);
    in.read((char*)p, x * sizeof(uint64_t));

//    std::cout << metadata_vec << std::endl;
    d_options.number_of_bins = metadata_vec[0];
    d_options.kmer_size = metadata_vec[2];

    return true;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
//
int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    DisOptions d_options;
    setupArgumentParser(parser, d_options);

    ArgumentParser::ParseResult res = parseCommandLine(d_options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (!read_filter_metadata(d_options))
        return 1;

    if (!verify_indices_dir(d_options.indices_dir,d_options.number_of_bins))
        return 1;

    if (!check_read_files(d_options))
        return 1;

    try
    {
        configure_d_mapper(d_options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
