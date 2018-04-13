// ==========================================================================
//                                 d_bloom_filter.h
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

#include <sdsl/bit_vectors.hpp>
#include <valarray>
#include <algorithm>

namespace seqan
{
    template<typename TString = Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna, SimpleShape> TShape;

        SeqAnBloomFilter(uint32_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBits(vec_size),
                        _filterVector(sdsl::bit_vector(vec_size, 0))

        {
            _init();
        }

        SeqAnBloomFilter(const char *fileName, uint32_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBits(vec_size)
        {
            _init();
            if (!sdsl::load_from_file(_filterVector, fileName))
            {
                std::cerr << "File \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }
            if (vec_size != _filterVector.bit_size())
            {
                std::cerr << "Size mismatch: \n\t Loaded file \t" << _filterVector.bit_size()
                << " bits\n\t expected\t" << vec_size << std::endl;
                exit(1);
            }
        }

        SeqAnBloomFilter(const char *fileName)
        {
            if (!sdsl::load_from_file(_filterVector, fileName))
            {
                std::cerr << "File \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }
            _noOfBits = _filterVector.bit_size();
            _getMetadata();
            _init();
        }

        void addKmers(TString const & text, uint32_t const & binNo)
        {
            _addKmers(text, binNo);
        }

        // ----------------------------------------------------------------------------
        // Function clearBins()
        // ----------------------------------------------------------------------------
        void clearBins(std::vector<uint32_t> & bins2clear, uint32_t & threadsCount)
        {
            std::vector<std::future<void>> tasks;

            uint64_t batchSize = _noOfHashPos/threadsCount;
            if(batchSize * threadsCount < _noOfHashPos) ++batchSize;

            for (uint32_t taskNo = 0; taskNo < threadsCount; ++taskNo)
            {
                tasks.emplace_back(std::async([=] {
                    for (uint64_t hashBlock=taskNo*batchSize; hashBlock < _noOfHashPos && hashBlock < (taskNo +1) * batchSize; ++hashBlock)
                    {
                        uint64_t vecPos = hashBlock * _blockBitSize;
                        for(uint32_t binNo : bins2clear)
                        {
                            _filterVector[vecPos + binNo] = false;
                        }
                    }
                }));
            }
            for (auto &&task : tasks)
            {
                task.get();
            }
        }
        // ----------------------------------------------------------------------------
        // Function addFastaFile()
        // ----------------------------------------------------------------------------
        void addFastaFile(CharString const & fastaFile, uint32_t const & binNo)
        {
            CharString id;
            IupacString seq;

            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(fastaFile)))
            {
                CharString msg = "Unable to open contigs File: ";
                append (msg, fastaFile);
                throw toCString(msg);
            }
            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                if(length(seq) < _kmerSize)
                    continue;
                addKmers(seq, binNo);
            }
            close(seqFileIn);
        }

        //save case sdsl
        bool save(const char *fileName)
        {
            _setMetadata();
            return sdsl::store_to_file(_filterVector, fileName);
        }

        double size_mb()
        {
            return sdsl::size_in_mega_bytes(_filterVector);
        }
        
        inline void whichBins(std::vector<bool> & selected, TString const & text, uint16_t const & threshold) const
        {
            uint16_t possible = length(text) - _kmerSize + 1;

            std::vector<uint16_t> counts(_noOfBins, 0);
            std::vector<uint64_t> kmerHashes(possible, 0);

            TShape kmerShape;
            resize(kmerShape, _kmerSize);
            hashInit(kmerShape, begin(text));
            auto it = begin(text);
            for (uint32_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = hashNext(kmerShape, it);
                ++it;
            }

            for (uint64_t kmerHash : kmerHashes)
            {
                std::vector<uint64_t> vecIndices =_preCalcValues;
                for(uint8_t i = 0; i < _noOfHashFunc ; i++)
                {
                    vecIndices[i] *= kmerHash;
                    getHashValue(vecIndices[i]);
                }
                uint32_t binNo = 0;
                for (uint16_t batchNo = 0; batchNo < _binIntWidth; ++batchNo)
                {
                    binNo = batchNo * INT_WIDTH;
                    uint64_t tmp = _filterVector.get_int(vecIndices[0], INT_WIDTH);
                    for(uint8_t i = 1; i < _noOfHashFunc;  i++)
                    {
                        tmp &= _filterVector.get_int(vecIndices[i], INT_WIDTH);
                    }

                    if (tmp ^ (1ULL<<(INT_WIDTH-1)))
                    {
                        while (tmp > 0)
                        {
                            uint64_t step = sdsl::bits::lo(tmp);
                            binNo += step;
                            ++step;
                            tmp >>= step;
                            ++counts[binNo];
                            ++binNo;
                        }
                    }
                    else
                    {
                        ++counts[binNo + INT_WIDTH - 1];
                    }
                    for(uint8_t i = 0; i < _noOfHashFunc ; i++)
                    {
                        vecIndices[i] += INT_WIDTH;
                    }
                }
            }

            for(uint32_t binNo=0; binNo < _noOfBins; ++binNo)
            {
                if(counts[binNo] >= threshold)
                    selected[binNo] = true;
            }
        }

        std::vector<bool> whichBins(TString const & text, uint16_t const & threshold) const
        {
            std::vector<bool> selected(_noOfBins, false);
            whichBins(selected, text, threshold);
            return selected;
        }

        uint32_t getNumberOfBins()
        {
            return _noOfBins;
        }

        uint8_t getKmerSize()
        {
            return _kmerSize;
        }

    private:

        void _init()
        {
            _binIntWidth = std::ceil((float)_noOfBins / INT_WIDTH);
            _blockBitSize = _binIntWidth * INT_WIDTH;
            _noOfHashPos = (_noOfBits - filterMetadataSize) / _blockBitSize;

            _preCalcValues.resize(_noOfHashFunc);
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
                _preCalcValues[i] = i ^  (_kmerSize * _seedValue);
        }
        void _getMetadata()
        {
            //-------------------------------------------------------------------
            //|              bf              | n_bins | n_hash_func | kmer_size |
            //-------------------------------------------------------------------
            uint64_t metadataStart = _noOfBits - filterMetadataSize;

            _noOfBins = _filterVector.get_int(metadataStart);
            _noOfHashFunc = _filterVector.get_int(metadataStart+64);
            _kmerSize = _filterVector.get_int(metadataStart+128);
        }

        void _setMetadata()
        {
            // -------------------------------------------------------------------
            // |              bf              | n_bins | n_hash_func | kmer_size |
            // -------------------------------------------------------------------
            uint64_t metadataStart = _noOfBits - filterMetadataSize;

            _filterVector.set_int(metadataStart, _noOfBins);
            _filterVector.set_int(metadataStart + 64, _noOfHashFunc);
            _filterVector.set_int(metadataStart + 128, _kmerSize);
        }

        template<typename TInt>
        bool _isBitSet(TInt num, uint8_t bit) const
        {
            return 1 == ( (num >> bit) & 1);
        }


        inline void getHashValue(uint64_t & vecIndex) const
        {
            vecIndex ^= vecIndex >> _shiftValue;
            vecIndex %= _noOfHashPos;
            vecIndex *= _blockBitSize;
        }

        void _insertKmer(uint64_t & kmerHash, uint32_t const & batchOffset)
        {
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                uint64_t vecIndex = _preCalcValues[i] * kmerHash;
                getHashValue(vecIndex);
                vecIndex += batchOffset;
                _filterVector[vecIndex] = 1;
            }
        }

        void _addKmers(TString const & text, uint32_t const & binNo)
        {
            TShape kmerShape;
            resize(kmerShape, _kmerSize);
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                _insertKmer(kmerHash, binNo);
            }
        }

        uint32_t                 _noOfBins;
        uint8_t                  _noOfHashFunc;
        uint8_t                  _kmerSize;
        uint16_t                 _binIntWidth;
        uint32_t                 _blockBitSize;

        //sizes in diferent units
        uint64_t                _noOfBits;
        uint64_t                _noOfHashPos;
        sdsl::bit_vector        _filterVector;

        std::vector<uint64_t>   _preCalcValues;
        uint64_t const          _shiftValue = 27;
        uint64_t const          _seedValue = 0x90b45d39fb6da1fa;
    };
}
