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
#include "d_kdx_filter.h"
#include "d_bloom_filter.h"
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
void setupArgumentParser(ArgumentParser & parser, DisOptions const & disOptions)
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

    // Setup output disOptions.
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
    setDefaultValue(parser, "read-group", disOptions.readGroup);
    addOption(parser, ArgParseOption("sm", "secondary-matches", "Specify whether to output secondary matches in \
                                     the XA tag of the primary alignment, as separate \
                                     secondary records, or to omit them.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "secondary-matches", disOptions.secondaryMatchesList);
    setDefaultValue(parser, "secondary-matches", disOptions.secondaryMatchesList[disOptions.secondaryMatches]);

    // Keep legacy option to display an error
    addOption(parser, ArgParseOption("sa", "secondary-alignments", "This option has been renamed to 'secondary-matches'.",
                                     ArgParseOption::STRING));
    hideOption(parser, "sa");
    addOption(parser, ArgParseOption("as", "align-secondary", "Compute and output co- and suboptimal \
                                     match alignments. Ignored if '-sa omit' is set."));

    addOption(parser, ArgParseOption("ra", "rabema-alignments", "Output alignments compatible with the \
                                     Read Alignment BEnchMArk (RABEMA)."));

    addOption(parser, ArgParseOption("sk", "skip-sam-headers", "Skip writing SQ: headers to SAM output (works only with SAM format)."));

    // Setup mapping disOptions.
    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("e", "error-rate", "Consider alignments within this percentual number of errors. \
                                     Increase this threshold to increase the number of mapped reads. \
                                     Decrease this threshold to decrease the runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", 100.0 * disOptions.errorRate);

    addOption(parser, ArgParseOption("s", "strata-rate", "Consider suboptimal alignments within this percentual number \
                                     of errors from the optimal alignment. Increase this threshold to increase \
                                     the number of alternative alignments at the expense of runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "strata-rate", "0");
    setMaxValue(parser, "strata-rate", "10");
    setDefaultValue(parser, "strata-rate", 100.0 * disOptions.strataRate);

    addOption(parser, ArgParseOption("y", "sensitivity", "Sensitivity with respect to edit distance. \
                                     Full sensitivity guarantees to find all considered alignments \
                                     but increases runtime, low sensitivity decreases runtime by \
                                     breaking such guarantee.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "sensitivity", disOptions.sensitivityList);
    setDefaultValue(parser, "sensitivity", disOptions.sensitivityList[disOptions.sensitivity]);

    // Setup paired-end mapping disOptions.
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
    setDefaultValue(parser, "indel-rate", 100.0 * disOptions.indelRate);

    addOption(parser, ArgParseOption("ni", "no-indels", "Turn off the rescue of unaligned ends containing indels."));

    //    addOption(parser, ArgParseOption("lo", "library-orientation", "Expected orientation of the segments in the library.",
    //                                     ArgParseOption::STRING));
    //    setValidValues(parser, "library-orientation", disOptions.libraryOrientationList);
    //    setDefaultValue(parser, "library-orientation", disOptions.libraryOrientationList[disOptions.libraryOrientation]);

    //    addOption(parser, ArgParseOption("la", "anchor", "Anchor one read and verify its mate."));

    // Setup performance disOptions.
    addSection(parser, "Performance Options");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", disOptions.threadsCount);

    addOption(parser, ArgParseOption("rb", "reads-batch", "Specify the number of reads to process in one batch.",
                                     ArgParseOption::INTEGER));
    
    setMinValue(parser, "reads-batch", "1000");
    setMaxValue(parser, "reads-batch", "5000000");
    setDefaultValue(parser, "reads-batch", disOptions.readsCount);
    hideOption(getOption(parser, "reads-batch"));

    // Setup Distributed mapper disOptions.
    addSection(parser, "Distributed mapper Options");
//    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
//                                     ArgParseOption::INTEGER));
//    setMinValue(parser, "number-of-bins", "1");
//    setMaxValue(parser, "number-of-bins", "1024");
//    setDefaultValue(parser, "number-of-bins", disOptions.numberOfBins);

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", disOptions.filterTypeList);
    setDefaultValue(parser, "filter-type",  disOptions.filterTypeList[disOptions.filterType]);


    addOption(parser, ArgParseOption("fi", "bloom-filter", "The path to a bloom filter. Default: will look for bloom.filter file inside the indices directory.", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "bloom-filter", "filter");
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(DisOptions & disOptions, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse indexed genome input file.
    getArgumentValue(disOptions.IndicesDirectory, parser, 0);

    // Append trailing slash if it doesn't exist.
    appendTrailingSlash(disOptions.IndicesDirectory);

    // Parse read input files.
    switch (getArgumentValueCount(parser, 1))
    {
        case 1:
            getArgumentValue(disOptions.readsFile.i1, parser, 1, 0);
            disOptions.singleEnd = true;
            break;
        case 2:
            getArgumentValue(disOptions.readsFile.i1, parser, 1, 0);
            getArgumentValue(disOptions.readsFile.i2, parser, 1, 1);
            disOptions.singleEnd = false;
            break;
        default:
            std::cerr << getAppName(parser) << ": Too many arguments!" << std::endl;
            return ArgumentParser::PARSE_ERROR;
    }

    // Parse output file.
    getOptionValue(disOptions.superOutputFile, parser, "output-file");

    // Parse output format.
    CharString outputFormat;
    if (getOptionValue(outputFormat, parser, "output-format"))
    {
        addLeadingDot(outputFormat);
        guessFormatFromFilename(outputFormat, disOptions.outputFormat);
    }
    else
        assign(disOptions.outputFormat, Sam());

#if SEQAN_HAS_ZLIB
    getOptionValue(disOptions.uncompressedBam, parser, "uncompressed-bam");
#endif

    // Parse output disOptions.
    getOptionValue(disOptions.readGroup, parser, "read-group");
    getOptionValue(disOptions.secondaryMatches, parser, "secondary-alignments", disOptions.secondaryMatchesList);
    getOptionValue(disOptions.rabema, parser, "rabema-alignments");

    getOptionValue(disOptions.alignSecondary, parser, "align-secondary");
    if (isSet(parser, "secondary-alignments"))
    {
        std::cerr << getAppName(parser) << ": The 'secondary-alignments' option has been renamed to 'secondary-matches'!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    if (disOptions.alignSecondary && disOptions.secondaryMatches == OMIT)
    {
        disOptions.alignSecondary = false;
        std::cerr << getAppName(parser) << ": WARNING, ignoring '-as' as '-sa omit' is set." << std::endl;
    }

    if (isSet(parser, "skip-sam-headers")) disOptions.skipSamHeader = true;

    // Parse mapping disOptions.
    unsigned errorRate;
    if (getOptionValue(errorRate, parser, "error-rate"))
        disOptions.errorRate = errorRate / 100.0;

    unsigned strataRate;
    if (getOptionValue(strataRate, parser, "strata-rate"))
        disOptions.strataRate = strataRate / 100.0;

    getOptionValue(disOptions.sensitivity, parser, "sensitivity", disOptions.sensitivityList);

    // Parse paired-end mapping disOptions.
    getOptionValue(disOptions.libraryLength, parser, "library-length");
    getOptionValue(disOptions.libraryDev, parser, "library-deviation");
    //    getOptionValue(disOptions.libraryOrientation, parser, "library-orientation", disOptions.libraryOrientationList);

    unsigned indelRate;
    if (getOptionValue(indelRate, parser, "indel-rate"))
        disOptions.indelRate = indelRate / 100.0;

    disOptions.verifyMatches = !isSet(parser, "no-indels");

    // Parse performance disOptions.
    getOptionValue(disOptions.threadsCount, parser, "threads");
    getOptionValue(disOptions.readsCount, parser, "reads-batch");

//    // Parse Distributed mapper options
//    getOptionValue(disOptions.numberOfBins, parser, "number-of-bins");
//    if (isSet(parser, "number-of-bins")) getOptionValue(disOptions.numberOfBins, parser, "number-of-bins");

    // Parse contigs index prefix.
    getOptionValue(disOptions.filterFile, parser, "bloom-filter");
    if (!isSet(parser, "bloom-filter"))
    {
        disOptions.filterFile = disOptions.IndicesDirectory;
        append(disOptions.filterFile, "bloom.filter");
    }

    getOptionValue(disOptions.filterType, parser, "filter-type", disOptions.filterTypeList);

    if (isSet(parser, "verbose")) disOptions.verbose = 1;
    if (isSet(parser, "very-verbose")) disOptions.verbose = 2;

    // Get version.
    disOptions.version = getVersion(parser);

    // Get command line.
    for (int i = 0; i < argc; i++)
    {
        append(disOptions.commandLine, argv[i]);
        appendValue(disOptions.commandLine, ' ');
    }
    eraseBack(disOptions.commandLine);

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureDisMapper(DisOptions & disOptions,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    if (disOptions.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawnDisMapper<TContigsSize, TContigsLen, uint32_t>(disOptions, threading, sequencing, distance);
    }
    else
    {
        spawnDisMapper<TContigsSize, TContigsLen, uint64_t>(disOptions, threading, sequencing, distance);
    }
}

template <typename TContigsSize, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureDisMapper(DisOptions & disOptions,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    if (disOptions.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureDisMapper<TContigsSize, uint32_t>(disOptions, threading, sequencing, distance);
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        configureDisMapper<TContigsSize, uint64_t>(disOptions, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureDisMapper(DisOptions & disOptions,
                        TThreading const & threading,
                        TSequencing const & sequencing,
                        TSeedsDistance const & distance)
{
    disOptions.contigsMaxLength = 0;
    disOptions.contigsSize = 0;
    disOptions.contigsSum = 0;
    disOptions.contigOffsets.resize(disOptions.numberOfBins, 0);
    // We aggregate individual limit here to configure the dis_mapper limits
    for (uint32_t i=0; i < disOptions.numberOfBins; ++i)
    {
        disOptions.contigOffsets[i] = disOptions.contigsSize;
        Options options = disOptions;
        appendFileName(options.contigsIndexFile, disOptions.IndicesDirectory, i);
        if (!openContigsLimits(options))
            throw RuntimeError("Error while opening contig limits file.");

       disOptions.contigsMaxLength   = std::max(options.contigsMaxLength, disOptions.contigsMaxLength);
       disOptions.contigsSize       += options.contigsSize;
       disOptions.contigsSum        += options.contigsSum;
    }

    if (disOptions.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureDisMapper<uint8_t>(disOptions, threading, sequencing, distance);
    }
    else if (disOptions.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureDisMapper<uint16_t>(disOptions, threading, sequencing, distance);
    }
    else
    {
        configureDisMapper<uint32_t>(disOptions, threading, sequencing, distance);
    }
}

template <typename TThreading, typename TSequencing>
void configureDisMapper(DisOptions & disOptions,
                        TThreading const & threading,
                        TSequencing const & sequencing)
{
    if (disOptions.sensitivity == FULL)
        return configureDisMapper(disOptions, threading, sequencing, EditDistance());
    else
        return configureDisMapper(disOptions, threading, sequencing, HammingDistance());
}

template <typename TThreading>
void configureDisMapper(DisOptions & disOptions,
                        TThreading const & threading)
{
    if (disOptions.singleEnd)
        configureDisMapper(disOptions, threading, SingleEnd());
    else
        configureDisMapper(disOptions, threading, PairedEnd());
}

void configureDisMapper(DisOptions & disOptions)
{
#ifdef _OPENMP
    if (disOptions.threadsCount > 1)
        configureDisMapper(disOptions, Parallel());
    else
#endif
        configureDisMapper(disOptions, Serial());
}

// ----------------------------------------------------------------------------
// Function checkReadFiles()
// ----------------------------------------------------------------------------
//
bool checkReadFiles(DisOptions const &  disOptions)
{
    // check if read file(s) exist(s)
    SeqFileIn seqFile;
    if (!open(seqFile, toCString(disOptions.readsFile.i1)))
    {
        std::cerr << "Unable to open read file "<< toCString(disOptions.readsFile.i1) <<"!\n";
        return false;
    }
    close(seqFile);
    if (!disOptions.singleEnd && !open(seqFile, toCString(disOptions.readsFile.i2)))
    {
        std::cerr << "Unable to open read file "<< toCString(disOptions.readsFile.i2) <<"!\n";
        return false;
    }
    close(seqFile);
    return true;
}

// ----------------------------------------------------------------------------
// Function readFilterMetadata()
// ----------------------------------------------------------------------------
//
bool readFilterMetadata(DisOptions &  disOptions)
{
    std::ifstream  in(toCString(disOptions.filterFile), std::ios::in | std::ios::binary);
    uint64_t x = filterMetadataSize/8; //bits -> bytes

    sdsl::int_vector<64>  metadataVec(filterMetadataSize/64, 0); //bits -> uint64_t
    in.seekg(-x, in.end); // seek from end of file

    uint64_t* p  = &(metadataVec[0]);
    in.read((char*)p, x * sizeof(uint64_t));

//    std::cout << metadataVec << std::endl;
    disOptions.numberOfBins = metadataVec[0];
    disOptions.kmerSize = metadataVec[2];

    return true;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
//
int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    DisOptions disOptions;
    setupArgumentParser(parser, disOptions);

    ArgumentParser::ParseResult res = parseCommandLine(disOptions, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (!readFilterMetadata(disOptions))
        return 1;

    if (!verifyIndicesDir(disOptions.IndicesDirectory, disOptions.numberOfBins))
        return 1;

    if (!checkReadFiles(disOptions))
        return 1;

    try
    {
        configureDisMapper(disOptions);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
