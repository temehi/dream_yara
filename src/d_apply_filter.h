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

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class FilterAppOptions
// --------------------------------------------------------------------------

struct FilterAppOptions
{
public:
    CharString              filter_file;
    Pair<CharString>        reads_file;
    std::string             output_file;
    float                   error_rate;
    unsigned                reads_count;
    unsigned                threads_count;

    double                  load_filter_time;
    double                  filter_reads_time;

    uint32_t                kmer_size;
    uint32_t                number_of_bins;
    uint64_t                passed_reads_count;

    FilterType              filter_type = BLOOM;

    bool stat_only;


    std::vector<std::string> filter_typeList = {"bloom", "kmer_direct", "none"};
    unsigned            verbose;


    FilterAppOptions() :
        error_rate(0.05f),
        reads_count(100000),
        threads_count(1),
        load_filter_time(0.0),
        filter_reads_time(0.0),
        kmer_size(20),
        number_of_bins(64),
        passed_reads_count(0),
        stat_only(false),
        verbose(0)
    {
#ifdef _OPENMP
        threads_count = std::thread::hardware_concurrency();
#endif
    }

    uint16_t get_threshold(uint16_t readLen)
    {
        uint16_t maxError = error_rate * readLen;

        // same as readLen - kmer_size + 1 - (maxError * kmer_size);
        if (kmer_size * (1 + maxError) > readLen)
            return 0;

        return readLen - kmer_size * (1 + maxError) + 1;
    }

};

// ----------------------------------------------------------------------------
// Class FilterApp
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct FilterApp
{
    // typedef typename If<IsSameType<TSequencing, PairedEnd>,
    //                     Pair<SeqFileIn>, SeqFileIn>::Type           TReadsFileIn;
    typedef typename TConfig::TThreading                         TThreading;
    typedef SeqStore<void, YaraReadsConfig>                      TReads;
    typedef PrefetchedFile<SeqFileIn, TReads, TThreading>        TReadsFile;

    FilterAppOptions &          options;
    Timer<double>                       timer;
    Stats<double>                       stats;

    TReads             reads;
    TReadsFile         reads_file;

    FilterApp(FilterAppOptions & options) :
        options(options),
        reads_file(options.reads_count)
    {};
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function write_filter_results()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TFilter>
inline void write_filter_results(FilterApp<TSpec, TConfig> & me, TFilter const & filter)
{
    start(me.timer);

    std::map<uint32_t, std::vector<bool>> filter_result;
    std::vector<bool> empty_result(me.options.number_of_bins, false);

    uint32_t number_of_reads = getReadsCount(me.reads.seqs);
    uint16_t avg_read_length = lengthSum(me.reads.seqs) / (number_of_reads * 2);
    uint16_t average_threshold = me.options.get_threshold(avg_read_length);


    if (average_threshold == 0)
    {
        std::cerr <<"[WARNING!] 0 k-mer is required to filter a read!\n";
        std::cerr <<"All reads will pass filteration and be mapped everywhere.\n ";
        std::cerr <<"This will be extremly slow.\n ";
        std::cerr <<"Choose an approprate error rate based on kmer size and read length\n";
    }

    uint32_t num_threads = me.options.threads_count;
    uint32_t batch_size = number_of_reads/num_threads;
    if(batch_size * num_threads < number_of_reads) ++batch_size;

    std::vector<std::future<void>> tasks;

    for (uint32_t read_id = 0; read_id < number_of_reads; ++read_id)
    {
        filter_result[read_id] = empty_result;
    }

    for (uint32_t task_num = 0; task_num < num_threads; ++task_num)
    {
        tasks.emplace_back(std::async([=, &me, &filter, &filter_result] {
            uint16_t threshold = 0;
            for (uint32_t read_id = task_num*batch_size; read_id < number_of_reads && read_id < (task_num +1) * batch_size; ++read_id)
            {
                threshold = me.options.get_threshold(length(me.reads.seqs[read_id]));
                select(filter, filter_result[read_id], me.reads.seqs[read_id], threshold);
                select(filter, filter_result[read_id], me.reads.seqs[read_id + number_of_reads], threshold);
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }
    stop(me.timer);
    me.options.filter_reads_time += getValue(me.timer);

    // get count data
    for (uint32_t read_id = 0; read_id < number_of_reads; ++read_id)
    {
        me.options.passed_reads_count +=  std::count(filter_result[read_id].begin(), filter_result[read_id].end(), true);
    }

    if (me.options.stat_only)
        return;
    
    std::cerr << "Writing results to a file ..." << std::endl;
    // get write the filter result
    std::ofstream filter_ostream(me.options.output_file, std::ios::binary);
    for (uint32_t read_id = 0; read_id < number_of_reads; ++read_id)
    {
        filter_ostream  << me.reads.names[read_id];
        for (uint32_t bin_num=0; bin_num<me.options.number_of_bins; ++bin_num)
        {
            filter_ostream  << "," << filter_result[read_id][bin_num];
        }
        filter_ostream  << "\n";
    }
    filter_ostream.close();
}

// ----------------------------------------------------------------------------
// Function load_reads()
// ----------------------------------------------------------------------------
// Loads one block of reads.

template <typename TSpec, typename TConfig>
inline void load_reads(FilterApp<TSpec, TConfig> & me)
{
    start(me.timer);

    readRecords(me.reads, me.reads_file);

    // Append reverse complemented reads.
    appendReverseComplement(me.reads);

    stop(me.timer);

    me.stats.loadReads += getValue(me.timer);
    me.stats.loadedReads += getReadsCount(me.reads.seqs);

    if (me.options.verbose > 1)
    {
        std::cerr << "Loading reads:\t\t\t" << me.timer << std::endl;
        std::cerr << "Reads count:\t\t\t" << getReadsCount(me.reads.seqs) << std::endl;
    }
}


// ----------------------------------------------------------------------------
// Function run_filter_app()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TFilter>
inline void run_filter_app(FilterApp<TSpec, TConfig> & me, TFilter const & filter)
{
    if (!open(me.reads_file, toCString(me.options.reads_file.i1)))
        throw RuntimeError("Error while opening reads file.");

    while (true)
    {
        if (me.options.verbose > 1) printRuler(std::cerr);
        load_reads(me);
        if (empty(me.reads.seqs)) break;

        write_filter_results(me, filter);
        clear(me.reads);
    }
    close(me.reads_file);
}
