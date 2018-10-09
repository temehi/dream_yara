// ==========================================================================
//                                 d_mapper.h
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


#ifndef APP_YARA_DIS_MAPPER_H_
#define APP_YARA_DIS_MAPPER_H_

#ifdef FILTER_SKIP_KMER
    const int _FILTER_SKIP_KMER = FILTER_SKIP_KMER;
#else
    const int _FILTER_SKIP_KMER = 0;
#endif

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class DisOptions
// --------------------------------------------------------------------------

struct DisOptions : public Options
{
public:
    CharString              indices_dir;
    CharString              filter_file;
    CharString              super_output_file;

    double                  load_filter      = 0.0;
    double                  filter_reads     = 0.0;
    double                  copy_reads       = 0.0;
    double                  copy_alignments  = 0.0;
    double                  move_cigars      = 0.0;

    bool                    skip_sam_header = false;

    uint32_t                kmer_size = 20;
    uint32_t                number_of_bins = 64;

    uint32_t                cur_bin_number = 0;
    uint64_t                filtered_reads = 0;
    std::vector<uint32_t>   contig_offsets;

    std::vector<std::vector<uint32_t>>          orig_read_id_map;
    std::map<uint32_t, String<CigarElement<>>>  collected_cigars;

    FilterType               filter_type = BLOOM;
    std::vector<std::string> filter_typeList = {"bloom", "kmer_direct", "none"};

    uint32_t get_contig_offsets()
    {
        return contig_offsets[cur_bin_number];
    }

    uint16_t get_threshold(uint16_t read_len)
    {
        uint16_t max_error = errorRate * read_len;

        // same as read_len - kmer_size + 1 - (max_error * kmer_size);
        if (kmer_size * (1 + max_error) > read_len)
            return 0;

        return (read_len - kmer_size * (1 + max_error) + 1)/_FILTER_SKIP_KMER;
    }

};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function append_stats()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void append_stats(Mapper<TSpec, TMainConfig> & d_mapper, Mapper<TSpec, TConfig> & child_mapper)
{
    d_mapper.stats.loadContigs    += child_mapper.stats.loadContigs;
    d_mapper.stats.loadReads      += child_mapper.stats.loadReads;
    d_mapper.stats.collectSeeds   += child_mapper.stats.collectSeeds;
    d_mapper.stats.findSeeds      += child_mapper.stats.findSeeds;
    d_mapper.stats.classifyReads  += child_mapper.stats.classifyReads;
    d_mapper.stats.rankSeeds      += child_mapper.stats.rankSeeds;
    d_mapper.stats.extendHits     += child_mapper.stats.extendHits;
    d_mapper.stats.sortMatches    += child_mapper.stats.sortMatches;
    d_mapper.stats.compactMatches += child_mapper.stats.compactMatches;
    d_mapper.stats.selectPairs    += child_mapper.stats.selectPairs;
    d_mapper.stats.verifyMatches  += child_mapper.stats.verifyMatches;
    d_mapper.stats.alignMatches   += child_mapper.stats.alignMatches;
    d_mapper.stats.writeMatches   += child_mapper.stats.writeMatches;
    d_mapper.stats.rescuedReads   += child_mapper.stats.rescuedReads;
}


// ----------------------------------------------------------------------------
// Function copy_matches()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copy_matches(Mapper<TSpec, TMainConfig> & d_mapper, Mapper<TSpec, TConfig> & child_mapper, DisOptions & d_options)
{
    start(d_mapper.timer);
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatch             TMatch;
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading         TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatchesAppender   TMatchesAppender;

    TMatchesAppender appender(d_mapper.matchesByCoord);

    uint32_t match_count = length(child_mapper.matchesByCoord);
    for (uint32_t i = 0; i < match_count; ++i)
    {
        TMatch currentMatch;
        uint32_t read_id         = child_mapper.matchesByCoord[i].readId;
        uint32_t orig_read_id     = read_id;
        if(d_options.filter_type != NONE)
            orig_read_id = d_options.orig_read_id_map[d_options.cur_bin_number][read_id];

        currentMatch.readId        = orig_read_id;
        currentMatch.contigId      = child_mapper.matchesByCoord[i].contigId + d_options.get_contig_offsets();
        currentMatch.isRev         = child_mapper.matchesByCoord[i].isRev;
        currentMatch.contigBegin   = child_mapper.matchesByCoord[i].contigBegin;
        currentMatch.contigEnd     = child_mapper.matchesByCoord[i].contigEnd;
        currentMatch.errors        = child_mapper.matchesByCoord[i].errors;
        appendValue(appender, currentMatch, Generous(), TThreading());
    }
    stop(d_mapper.timer);
    d_options.copy_alignments += getValue(d_mapper.timer);
}

// ----------------------------------------------------------------------------
// Function display Cigars()
// ----------------------------------------------------------------------------
inline CharString cigars2String(String<CigarElement<> > const & cigar_str)
{
    CharString cs;
    for (CigarElement<> const & cigar_element : cigar_str)
    {
        appendNumber(cs, cigar_element.count);
        appendValue(cs, cigar_element.operation);
    }
    return cs;
}

inline StringSet<CharString> cigarsSet2String(StringSet<String<CigarElement<> >, Segment< String<CigarElement<> >> > & cigar_str_set)
{
    StringSet<CharString> css;
    for (String<CigarElement<>> const & cigar_str : cigar_str_set)
    {
        appendValue(css, cigars2String(cigar_str));
    }
    return css;
}


// ----------------------------------------------------------------------------
// Function copy_cigars()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copy_cigars(Mapper<TSpec, TMainConfig> & d_mapper, Mapper<TSpec, TConfig> & child_mapper, DisOptions & d_options)
{
    start(d_mapper.timer);
    typedef typename MapperTraits<TSpec, TConfig>::TMatch             TMatch;
    uint32_t match_count = length(child_mapper.primaryMatches);
    for (uint32_t i = 0; i < match_count; ++i)
    {
        TMatch currentMatch = child_mapper.primaryMatches[i];
        uint32_t read_id         = currentMatch.readId;
        uint32_t orig_read_id     = read_id;
        if(d_options.filter_type != NONE)
            orig_read_id = d_options.orig_read_id_map[d_options.cur_bin_number][read_id];

        setMapped(d_mapper.ctx, orig_read_id);
        setMinErrors(d_mapper.ctx, orig_read_id, currentMatch.errors);

        if (!isPaired(d_mapper.ctx, orig_read_id) && isPaired(child_mapper.ctx, read_id)) // First paired match
        {
            setPaired(d_mapper.ctx, orig_read_id);
            d_mapper.primaryMatchesProbs[orig_read_id] = child_mapper.primaryMatchesProbs[read_id];
        }

        if(getMinErrors(d_mapper.ctx, orig_read_id) == currentMatch.errors)
        {
            d_options.collected_cigars[orig_read_id] = child_mapper.cigars[read_id];
        }
    }
    stop(d_mapper.timer);
    d_options.move_cigars += getValue(d_mapper.timer);
}

// ----------------------------------------------------------------------------
// Function transfer_cigars()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig>
inline void transfer_cigars(Mapper<TSpec, TMainConfig> & d_mapper, DisOptions & d_options)
{
    start(d_mapper.timer);

    typedef MapperTraits<TSpec, TMainConfig>         TTraits;
    typedef typename TTraits::TCigarsPos             TCigarsPos;
    typedef typename TTraits::TThreading             TThreading;

    resize(d_mapper.cigars.limits, getReadsCount(d_mapper.reads.seqs)+1, 0);
    for(auto iter = d_options.collected_cigars.begin(); iter != d_options.collected_cigars.end(); ++iter)
    {
        d_mapper.cigars.limits[iter->first + 1] = length(iter->second);
        append(d_mapper.cigarString, iter->second);
    }
    partialSum(d_mapper.cigars.limits, TThreading());
    assign(d_mapper.cigars.positions, prefix(d_mapper.cigars.limits, length(d_mapper.cigars.limits) - 1));

    // If only the primary matches were aligned, we use the identity modifier
    setHost(d_mapper.primaryCigars, d_mapper.cigars);
    setCargo(d_mapper.primaryCigars, d_mapper.primaryCigarPositions);
    assign(d_mapper.primaryCigarPositions, seqan::Range<TCigarsPos>(0, length(d_mapper.primaryMatches)), Exact());

    stop(d_mapper.timer);
    d_options.move_cigars += getValue(d_mapper.timer);
}


// ----------------------------------------------------------------------------
// Function _map_reads_impl()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig, typename TReadSeqs>
inline void _map_reads_impl(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & d_mapper, TReadSeqs & read_seqs, DisOptions & d_options)
{
    initReadsContext(me, read_seqs);
    initSeeds(me, read_seqs);

    collectSeeds<0>(me, read_seqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, read_seqs);
    collectSeeds<2>(me, read_seqs);
    findSeeds<0>(me, 1);
    findSeeds<0>(me, 2);
    rankSeeds(me);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<0>(me, 1);
    extendHits<0>(me, 2);
    clearSeeds(me);
    clearHits(me);

    initSeeds(me, read_seqs);
    collectSeeds<1>(me, read_seqs);
    findSeeds<1>(me, 1);
    collectSeeds<2>(me, read_seqs);
    findSeeds<1>(me, 2);
    rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 1.
    extendHits<1>(me, 1);
    extendHits<1>(me, 2);
    clearSeeds(me);
    clearHits(me);

    if (me.options.sensitivity > LOW)
    {
        initSeeds(me, read_seqs);
        collectSeeds<2>(me, read_seqs);
        findSeeds<2>(me, 2);
        rankSeeds(me);
        // TODO(esiragusa): filter out hits with distance < 2.
        extendHits<2>(me, 2);
        clearHits(me);
        clearSeeds(me);
    }
    aggregateMatches(me, read_seqs);
    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    copy_matches(d_mapper, me, d_options);
    copy_cigars(d_mapper, me, d_options);
    append_stats(d_mapper, me);
}

// ----------------------------------------------------------------------------
// Function clasify_loaded_reads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void clasify_loaded_reads(Mapper<TSpec, TMainConfig>  & d_mapper, TFilter const & filter, DisOptions & d_options)
{
    start(d_mapper.timer);

    uint32_t number_of_reads = getReadsCount( d_mapper.reads.seqs);
    uint16_t avgReadLen = lengthSum(d_mapper.reads.seqs) / (number_of_reads * 2);
    uint16_t threshold = d_options.get_threshold(avgReadLen);

    d_options.orig_read_id_map.clear();
    d_options.orig_read_id_map.resize(d_options.number_of_bins);

    if (threshold == 0)
    {
        std::cerr <<"[WARNING!] 0 k-mer is required to filter a read!\n";
        std::cerr <<"All reads will pass filteration and be mapped everywhere.\n ";
        std::cerr <<"This will be extremly slow.\n ";
        std::cerr <<"Choose an approprate error rate based on kmer size and read length\n";
        for (uint32_t readID = 0; readID < number_of_reads; ++readID)
        {
            for (uint32_t binNo = 0; binNo < d_options.number_of_bins; ++binNo)
            {
               d_options.orig_read_id_map[binNo].push_back(readID);
            }
        }
    }
    // if paired classify only one pair
    if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        number_of_reads = getPairsCount( d_mapper.reads.seqs);

    uint32_t numThr = d_options.threadsCount;
    uint32_t batchSize = number_of_reads/numThr;
    if(batchSize * numThr < number_of_reads) ++batchSize;

    std::vector<std::future<void>> tasks;


    for (uint32_t taskNo = 0; taskNo < numThr; ++taskNo)
    {
        tasks.emplace_back(std::async([=, &d_mapper, &d_options, &filter] {
            for (uint32_t readID = taskNo*batchSize; readID < number_of_reads && readID < (taskNo +1) * batchSize; ++readID)
            {
                std::vector<bool> selectedBins(d_options.number_of_bins, false);
                // select(filter, selectedBins, d_mapper.reads.seqs[readID], threshold);
                // select(filter, selectedBins, d_mapper.reads.seqs[readID + number_of_reads], threshold);
                select<Offset<_FILTER_SKIP_KMER>>(filter, selectedBins, d_mapper.reads.seqs[readID], threshold);
                select<Offset<_FILTER_SKIP_KMER>>(filter, selectedBins, d_mapper.reads.seqs[readID + number_of_reads], threshold);

                if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
                {
                    // select(filter, selectedBins, d_mapper.reads.seqs[readID + 2*number_of_reads], threshold);
                    // select(filter, selectedBins, d_mapper.reads.seqs[readID + 3*number_of_reads], threshold);
                    select<Offset<_FILTER_SKIP_KMER>>(filter, selectedBins, d_mapper.reads.seqs[readID + 2*number_of_reads], threshold);
                    select<Offset<_FILTER_SKIP_KMER>>(filter, selectedBins, d_mapper.reads.seqs[readID + 3*number_of_reads], threshold);
                }

                for (uint32_t binNo = 0; binNo < d_options.number_of_bins; ++binNo)
                {
                    if(selectedBins[binNo])
                    {
                        mtx.lock();
                        d_options.orig_read_id_map[binNo].push_back(readID);
                        mtx.unlock();
                    }
                }
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }

    stop(d_mapper.timer);
    d_options.filter_reads += getValue(d_mapper.timer);

    if (d_options.verbose > 1)
    {
        for (uint32_t binNo = 0; binNo < d_options.number_of_bins; ++binNo)
        {
            std::cerr << "bin " << binNo << "\t" << d_options.orig_read_id_map[binNo].size() << std::endl;
        }
    }
}


// ----------------------------------------------------------------------------
// Function load_filtered_reads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void load_filtered_reads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & d_mapper, DisOptions & d_options)
{

    start(d_mapper.timer);

    if(d_options.filter_type == NONE)
    {
        for (uint32_t i = 0; i< getReadSeqsCount(d_mapper.reads.seqs); ++i)
        {
            appendValue(me.reads.seqs, d_mapper.reads.seqs[i]);
        }
    }
    else
    {
        uint32_t number_of_reads = getReadsCount( d_mapper.reads.seqs);
        uint32_t number_of_iltered_reads = d_options.orig_read_id_map[d_options.cur_bin_number].size();

        //load forward reads
        for (uint32_t i = 0; i< number_of_iltered_reads; ++i)
        {
            uint32_t orgId = d_options.orig_read_id_map[d_options.cur_bin_number][i];
            appendValue(me.reads.seqs, d_mapper.reads.seqs[orgId]);
        }

        // if paired classify only one pair
        if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        {
            uint32_t numPairs = getPairsCount( d_mapper.reads.seqs);
            //load mates
            for (uint32_t i = 0; i< number_of_iltered_reads; ++i)
            {
                uint32_t orgId = d_options.orig_read_id_map[d_options.cur_bin_number][i];
                appendValue(me.reads.seqs, d_mapper.reads.seqs[orgId + numPairs]);
                d_options.orig_read_id_map[d_options.cur_bin_number].push_back(orgId + numPairs);
            }
            number_of_iltered_reads *= 2; //now we have twice the reads
        }

        //load reverse reads
        for (uint32_t i = 0; i< number_of_iltered_reads; ++i)
        {
            uint32_t orgId = d_options.orig_read_id_map[d_options.cur_bin_number][i];
            appendValue(me.reads.seqs, d_mapper.reads.seqs[orgId + number_of_reads]);
            d_options.orig_read_id_map[d_options.cur_bin_number].push_back(orgId + number_of_reads);
        }
    }
    stop(d_mapper.timer);
    d_options.copy_reads += getValue(d_mapper.timer);
    d_options.filtered_reads += getReadsCount(me.reads.seqs);
}

// ----------------------------------------------------------------------------
// Function map_reads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void map_reads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & d_mapper, DisOptions & d_options)
{
    _map_reads_impl(me, d_mapper, me.reads.seqs, d_options);
}

template <typename TSpec, typename TConfig, typename TMainConfig>
inline void run_d_mapper(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig> & d_mapper, DisOptions & d_options)
{
    load_filtered_reads(me, d_mapper, d_options);
    if (empty(me.reads.seqs)) return;
    loadContigs(me);
    loadContigsIndex(me);
    map_reads(me, d_mapper, d_options);
}


// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TContigsSum,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void spawnMapper(Options const & options,
                 Mapper<TSpec, TMainConfig> & d_mapper,
                 DisOptions & d_options,
                 TThreading const & /*threading*/,
                 TSequencing const & /*sequencing*/,
                 TSeedsDistance const & /*distance*/)
{

    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TConfig;
    Mapper<void, TConfig> mapper(options);
    run_d_mapper(mapper, d_mapper, d_options);
}
// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & d_mapper,
                     DisOptions & d_options,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawnMapper<TContigsSize, TContigsLen, uint32_t>(options, d_mapper, d_options, threading, sequencing, distance);
    }
    else
    {
        spawnMapper<TContigsSize, TContigsLen, uint64_t>(options, d_mapper, d_options, threading, sequencing, distance);
    }
}

template <typename TContigsSize,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & d_mapper,
                     DisOptions & d_options,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureMapper<TContigsSize, uint32_t>(options, d_mapper, d_options, threading, sequencing, distance);
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMapper<TContigsSize, uint64_t>(options, d_mapper, d_options, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TSpec, typename TMainConfig>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & d_mapper,
                     DisOptions & d_options)
{
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading       TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSequencing      TSequencing;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSeedsDistance   TSeedsDistance;

    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureMapper<uint8_t>(options, d_mapper, d_options, TThreading(), TSequencing(), TSeedsDistance());
    }
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureMapper<uint16_t>(options, d_mapper, d_options, TThreading(), TSequencing(), TSeedsDistance());
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMapper<uint32_t>(options, d_mapper, d_options, TThreading(), TSequencing(), TSeedsDistance());
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}


// ----------------------------------------------------------------------------
// Function loadAllContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadAllContigs(Mapper<TSpec, TConfig> & d_mapper, DisOptions & d_options)
{
    typedef typename MapperTraits<TSpec, TConfig>::TContigs          TContigs;

    start(d_mapper.timer);
    try
    {
        for (uint32_t i=0; i < d_options.number_of_bins; ++i)
        {
            TContigs tmp_contigs;
            CharString file_name;
            append_file_name(file_name, d_options.indices_dir, i);

            if (!open(tmp_contigs, toCString(file_name), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");
            append(d_mapper.contigs.seqs, tmp_contigs.seqs);
            append(d_mapper.contigs.names, tmp_contigs.names);
        }
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(d_mapper.timer);
    d_mapper.stats.loadContigs += getValue(d_mapper.timer);

    if (d_mapper.options.verbose > 1)
        std::cerr << "Loading reference:\t\t\t" << d_mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function rank_matches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void rank_matches(Mapper<TSpec, TConfig> & me, TReadSeqs const & read_seqs)
{
    typedef MapperTraits<TSpec, TConfig>                    TTraits;
    typedef typename TTraits::TMatch                        TMatch;
    typedef typename TTraits::TMatchesSet                   TMatchesSet;
    typedef typename TTraits::TMatchesViewSet               TMatchesViewSet;
    typedef typename Value<TMatchesSet const>::Type         TMatchesSetValue;
    typedef typename Value<TMatchesViewSet const>::Type     TMatchesViewSetValue;
    typedef typename Iterator<TMatchesViewSet const, Standard>::Type TMatchesViewSetIt;
    typedef typename Size<TReadSeqs>::Type                  TReadId;
    typedef typename Size<TMatchesSetValue>::Type           TMatchesSize;
    typedef std::uniform_int_distribution<TMatchesSize>     TMatchesRnd;

    start(me.timer);
    // Create a position modifier of the matches from the identity permutation.
    assign(me.matchesPositions, seqan::Range<TMatchesSize>(0, length(me.matchesByCoord)), Exact());
    setHost(me.matchesByErrors, me.matchesByCoord);
    setCargo(me.matchesByErrors, me.matchesPositions);

    // Bucket matches in the position modifier.
    setHost(me.matchesSetByErrors, me.matchesByErrors);
    assign(stringSetLimits(me.matchesSetByErrors), stringSetLimits(me.matchesSetByCoord), Exact());
    assign(stringSetPositions(me.matchesSetByErrors), stringSetPositions(me.matchesSetByCoord), Exact());

    // Sort matches by pairing info. iff possible

    // Sort matches by errors.
    forEach(me.matchesSetByErrors, sortMatches<TMatchesViewSetValue, Errors>, typename TTraits::TThreading());

    // Select all co-optimal matches.
    assign(me.optimalMatchesSet, me.matchesSetByErrors);
    clipMatches(me.optimalMatchesSet, countMatchesInBestStratum<TMatchesViewSetValue>, typename TTraits::TThreading());

    // Select all sub-optimal matches.
    assign(me.matchesSet, me.matchesSetByErrors);
    clipMatches(me.matchesSet, [&](TMatchesViewSetValue const & matches)
                {
                    if (empty(matches)) return TMatchesSize(0);

                    TReadId read_id = getMember(front(matches), ReadId());

                    return countMatchesInStrata(matches, getReadStrata<TMatch>(me.options, length(read_seqs[read_id])));
                },
                typename TTraits::TThreading());

    // Append an invalid match to matches by coord.
    resize(me.matchesByCoord, length(me.matchesByCoord) + 1, Exact());
    setInvalid(back(me.matchesByCoord));
    // Update matches by errors.
    resize(me.matchesPositions, length(me.matchesPositions) + 1, Exact());
    setPosition(me.matchesByErrors, length(me.matchesByErrors) - 1, length(me.matchesByCoord) - 1);

    // Initialize primary matches.
    setHost(me.primaryMatches, me.matchesByErrors);
    assign(me.primaryMatchesPositions, stringSetPositions(me.matchesSetByErrors), Exact());
    setCargo(me.primaryMatches, me.primaryMatchesPositions);

    // Choose primary matches among best matches.
    iterate(me.optimalMatchesSet, [&](TMatchesViewSetIt const & matchesIt)
            {
                // Use one generator per thread.
                std::default_random_engine generator;

                TReadId read_id = position(matchesIt, me.optimalMatchesSet);
                TMatchesViewSetValue const & matches = value(matchesIt);

                // Set unmapped reads as invalid.
                if (empty(matches))
                {
                    setPosition(me.primaryMatches, read_id, length(me.matchesByErrors) - 1);
                }
                // Choose match at random.
                else
                {
                    TMatchesRnd rnd(0, length(matches) - 1);
                    setPosition(me.primaryMatches, read_id, position(me.primaryMatches, read_id) + rnd(generator));
                }
            },
            Standard(), typename TTraits::TThreading());

    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);
    if (me.options.verbose > 1)
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Update mapped reads.
    transform(me.ctx.mapped, me.primaryMatches, isValid<typename TTraits::TMatchSpec>, typename TTraits::TThreading());

    if (me.options.verbose > 0)
    {
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;

        if (me.options.verbose > 1)
            std::cerr << "Mapped reads:\t\t\t" << mappedReads << std::endl;
    }

    if (IsSameType<typename TConfig::TSequencing, SingleEnd>::VALUE) return;

    // Update paired reads.
    if (me.options.verbose > 0)
    {
        unsigned long pairedReads = count(me.ctx.paired, true, typename TTraits::TThreading());
        me.stats.pairedReads += pairedReads;

        if (me.options.verbose > 1)
        {
            std::cerr << "Pairing time:\t\t\t" << me.timer << std::endl;
            std::cerr << "Paired reads:\t\t\t" << pairedReads << std::endl;
        }
    }
}

// ----------------------------------------------------------------------------
// Function open_output_file()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig>
inline void open_output_file(Mapper<TSpec, TConfig> & d_mapper, DisOptions & d_options)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef typename TTraits::TContigs              TContigs;
    typedef typename TTraits::TContigSeqs           TContigSeqs;
    typedef typename Value<TContigSeqs>::Type       TContigSeq;

    String<uint32_t> allContigLengths;

    start(d_mapper.timer);
    try
    {
        for (uint32_t i=0; i < d_options.number_of_bins; ++i)
        {
            TContigs tmp_contigs;
            String<uint32_t> tmp_contig_len;
            CharString file_name;
            append_file_name(file_name, d_options.indices_dir, i);

            if (!open(tmp_contigs, toCString(file_name), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");

            resize(tmp_contig_len, length(tmp_contigs.seqs));
            transform(tmp_contig_len, tmp_contigs.seqs, [](TContigSeq const & seq) { return length(seq); });
            append(allContigLengths, tmp_contig_len);

            append(d_mapper.contigs.names, tmp_contigs.names);
        }
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(d_mapper.timer);
    d_mapper.stats.loadContigs += getValue(d_mapper.timer);

    if (d_mapper.options.verbose > 1)
        std::cerr << "Loading reference:\t\t\t" << d_mapper.timer << std::endl;

    bool opened = false;

    if (empty(d_mapper.options.outputFile))
    {
        // Output to cout.
        if (d_mapper.options.uncompressedBam)
        {
            // Turn off BAM compression.
            setFormat(d_mapper.outputFile, d_mapper.options.outputFormat);
            opened = _open(d_mapper.outputFile, std::cout, Nothing(), False());
        }
        else
        {
            opened = open(d_mapper.outputFile, std::cout, d_mapper.options.outputFormat);
        }
    }
    else
    {
        // Output to file.
        opened = open(d_mapper.outputFile, toCString(d_mapper.options.outputFile), OPEN_WRONLY | OPEN_CREATE);
    }

    if (!opened) throw RuntimeError("Error while opening output file.");

    setContigNames(context(d_mapper.outputFile), d_mapper.contigs.names);

    // Fill contig lengths.
    resize(contigLengths(context(d_mapper.outputFile)), length(allContigLengths));
    assign(contigLengths(context(d_mapper.outputFile)), allContigLengths);

    typedef FileFormat<BamFileOut>::Type    TOutputFormat;
    TOutputFormat of;
    assign(of, Bam());

    if(d_mapper.outputFile.format.tagId == of.tagId || !d_options.skip_sam_header)
    {
        // Write header.
        BamHeader header;
        fillHeader(header, d_mapper.options);
        writeHeader(d_mapper.outputFile, header);
    }
}

// ----------------------------------------------------------------------------
// Function prepair_d_mapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void prepair_d_mapper(Mapper<TSpec, TMainConfig> & d_mapper, TFilter const & filter, DisOptions & d_options)
{
    initReadsContext(d_mapper, d_mapper.reads.seqs);
    setHost(d_mapper.cigars, d_mapper.cigarString);
    d_options.collected_cigars.clear();
    if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        resize(d_mapper.primaryMatchesProbs, getReadsCount(d_mapper.reads.seqs), 0.0, Exact());
    if(d_options.filter_type != NONE)
        clasify_loaded_reads(d_mapper, filter, d_options);
}

// ----------------------------------------------------------------------------
// Function finalize_d_mapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig>
inline void finalize_d_mapper(Mapper<TSpec, TMainConfig> & d_mapper, DisOptions & d_options)
{
    aggregateMatches(d_mapper, d_mapper.reads.seqs);
    rank_matches(d_mapper, d_mapper.reads.seqs);
    transfer_cigars(d_mapper, d_options);

    writeMatches(d_mapper);
    clearMatches(d_mapper);
    clearAlignments(d_mapper);
    clearReads(d_mapper);
}

// ----------------------------------------------------------------------------
// Function sorted_bins()
// ----------------------------------------------------------------------------
std::vector<uint32_t> sorted_bins(DisOptions const & d_options)
{

    std::vector<uint32_t> sorted_bin_index(d_options.number_of_bins);
    iota(sorted_bin_index.begin(), sorted_bin_index.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(sorted_bin_index.begin(), sorted_bin_index.end(),
         [&d_options](size_t i1, size_t i2) {return d_options.orig_read_id_map[i1].size() > d_options.orig_read_id_map[i2].size();});

    return sorted_bin_index;
}

// ----------------------------------------------------------------------------
// Function run_d_mapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void run_d_mapper(Mapper<TSpec, TMainConfig> & d_mapper, TFilter const & filter, DisOptions & d_options)
{
    configureThreads(d_mapper);

    // Open output file and write header.
    open_output_file(d_mapper, d_options);
    openReads(d_mapper);

    while (true)
    {
        if (d_mapper.options.verbose > 1) printRuler(std::cerr);
        loadReads(d_mapper);
        if (empty(d_mapper.reads.seqs)) break;

        prepair_d_mapper(d_mapper, filter, d_options);

        for (auto i: sorted_bins(d_options))
        {
            d_options.cur_bin_number = i;
            Options options = d_mapper.options;
            append_file_name(options.contigsIndexFile, d_options.indices_dir, i);
            if (!openContigsLimits(options))
                throw RuntimeError("Error while opening reference file.");
            configureMapper<TSpec, TMainConfig>(options, d_mapper, d_options);
        }

        finalize_d_mapper(d_mapper, d_options);
    }
    closeReads(d_mapper);
    closeOutputFile(d_mapper);
}

// ----------------------------------------------------------------------------
// Function spawn_d_mapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
typename TThreading, typename TSequencing, typename TSeedsDistance>
inline void spawn_d_mapper(DisOptions & d_options,
                           TThreading const & /* tag */,
                           TSequencing const & /* tag */,
                           TSeedsDistance const & /* tag */)
{
    d_options.outputFile = d_options.super_output_file;
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TMainConfig;
    Mapper<void, TMainConfig> d_mapper(d_options);

    Timer<double> timer;

    start(timer);
    start(d_mapper.timer);

    if (d_options.filter_type == BLOOM)
    {
        typedef BinningDirectory<InterleavedBloomFilter,
                                 BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;

        BinningDirectoriesIBF filter;
        retrieve(filter, toCString(d_options.filter_file));

        d_options.kmer_size = filter.kmerSize;

        stop(d_mapper.timer);
        d_options.load_filter += getValue(d_mapper.timer);
        run_d_mapper(d_mapper, filter, d_options);
    }
    else if (d_options.filter_type == KMER_DIRECT)
    {
        typedef BinningDirectory<DirectAddressing,
                                 BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesDA;
        BinningDirectoriesDA filter;
        retrieve(filter, toCString(d_options.filter_file));

        d_options.kmer_size = filter.kmerSize;

        stop(d_mapper.timer);
        d_options.load_filter += getValue(d_mapper.timer);
        run_d_mapper(d_mapper, filter, d_options);
    }
    else
    {
        // dummy filter in case of nofilter option
        typedef BinningDirectory<InterleavedBloomFilter,
                                 BDConfig<Dna, Normal, Uncompressed> > BinningDirectoriesIBF;

        BinningDirectoriesIBF filter(64, 3, 20, 1);

        stop(d_mapper.timer);
        d_options.load_filter += getValue(d_mapper.timer);
        run_d_mapper(d_mapper, filter, d_options);
    }
    stop(timer);
    if (d_mapper.options.verbose > 0)
    {
        double total = getValue(timer) / 100.0;

        std::cerr << "\nFilter loading time:\t\t" << d_options.load_filter << " sec" << "\t\t" << d_options.load_filter / total << " %" << std::endl;
        std::cerr << "Reads filtering time:\t\t" << d_options.filter_reads << " sec" << "\t\t" << d_options.filter_reads / total << " %" << std::endl;
        std::cerr << "Reads copying time:\t\t" << d_options.copy_reads << " sec" << "\t\t" << d_options.copy_reads / total << " %" << std::endl;
        std::cerr << "Alignments copying time:\t" << d_options.copy_alignments << " sec" << "\t\t" << d_options.copy_alignments / total << " %" << std::endl;
        std::cerr << "Cigars moving time:\t\t" << d_options.move_cigars << " sec" << "\t\t" << d_options.move_cigars / total << " %" << std::endl;

        printStats(d_mapper, timer);
		std::cerr << "Avg reads per bin:\t\t" << (double)d_options.filtered_reads / d_options.number_of_bins << std::endl;
	}
}

#endif  // #ifndef APP_YARA_MAPPER_H_
