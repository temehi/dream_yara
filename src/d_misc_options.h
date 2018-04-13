// ==========================================================================
//                                 d_misc_options.h
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
#include <cstdio>

#ifndef APP_YARA_MISC_OPTIONS_DIS_H_
#define APP_YARA_MISC_OPTIONS_DIS_H_
std::mutex mtx;

static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;


enum FilterType
{
    BLOOM, KMER_DIRECT, NONE
};

using namespace seqan;
class Semaphore
{
    std::mutex m;
    std::condition_variable cv;
    int count;

public:
    Semaphore(int n) : count{n} {}
    void notify()
    {
        std::unique_lock<std::mutex> l(m);
        ++count;
        cv.notify_one();
    }
    void wait()
    {
        std::unique_lock<std::mutex> l(m);
        cv.wait(l, [this]{ return count!=0; });
        --count;
    }
};

class Critical_section
{
    Semaphore &s;
public:
    Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
    ~Critical_section() { s.notify(); }
};

// ============================================================================
// Functions
// ============================================================================
// ----------------------------------------------------------------------------
// Function appendFileName()
// ----------------------------------------------------------------------------
inline void appendFileName(CharString & target, CharString const & source, uint32_t const i)
{
    target = source;
    append(target, std::to_string(i));
}

inline void appendFileName(CharString & target, uint32_t const i)
{
    CharString source = target;
    appendFileName(target, source, i);
}

// ----------------------------------------------------------------------------
// Function getExtensionWithLeadingDot()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString const>::Type
getExtensionWithLeadingDot(TString const & string)
{
    return suffix(string, firstOf(string, IsDot()));
}

// ----------------------------------------------------------------------------
// Function getFilesInDir()
// ----------------------------------------------------------------------------
inline void getFilesInDir(StringSet<CharString> & fileNames, CharString const directoryPath)
{
    DIR *dir;
    struct dirent *ent;
    struct stat st;

    dir = opendir(toCString(directoryPath));
    while ((ent = readdir(dir)) != NULL)
    {
        CharString fileName = ent->d_name;
        CharString fullFileName = directoryPath;
        append(fullFileName, "/");
        append(fullFileName, fileName);


        bool invalidFile = (fileName[0] == '.')
        || (stat(toCString(fullFileName), &st) == -1)
        || ((st.st_mode & S_IFDIR) != 0);

        if (!invalidFile)
            appendValue(fileNames, fileName);
    }
}

// ----------------------------------------------------------------------------
// Function getValidFilesInDir()
// ----------------------------------------------------------------------------
inline void getValidFilesInDir(StringSet<CharString> & fileNames,
                               CharString const directoryPath,
                               std::vector<std::string> const & validExtensions)
{
    StringSet<CharString>  allFileNames;
    getFilesInDir(allFileNames, directoryPath);
    for (uint32_t i = 0; i < length(allFileNames); ++i)
    {
        CharString ext = getExtensionWithLeadingDot(allFileNames[i]);

        auto it = std::find(validExtensions.begin(), validExtensions.end(), std::string(toCString(ext)));
        if (it != validExtensions.end())
            appendValue(fileNames, allFileNames[i]);
    }

}

// ----------------------------------------------------------------------------
// Function verifyIndicesDir()
// ----------------------------------------------------------------------------
inline bool verifyIndicesDir(CharString const directoryPath, uint32_t const numberOfBins)
{
    for (uint32_t i=0; i < numberOfBins; ++i)
    {
        CharString contigsLimitFile;
        appendFileName(contigsLimitFile, directoryPath, i);
        append(contigsLimitFile, ".txt.size");

        String<uint64_t> limits;

        if (!open(limits, toCString(contigsLimitFile), OPEN_RDONLY|OPEN_QUIET))
        {
            std::cerr << "No index for bin " << i << '\n';
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function verifyFnaFile()
// ----------------------------------------------------------------------------
inline bool verifyFnaFile(CharString const & fastaFile)
{
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(fastaFile)))
    {
        std::cerr << "Fasta file: " << fastaFile << " can not be found!\n" ;
        return false;
    }
    close(seqFileIn);
    return true;
}

// ----------------------------------------------------------------------------
// Function commonExtension()
// ----------------------------------------------------------------------------
inline std::string commonExtension(CharString const directoryPath, uint32_t const numberOfBins)
{
    std::vector<std::string> extensions =  SeqFileIn::getFileExtensions();

    uint32_t count = 0;
    for (auto ext : extensions)
    {
        count = 0;
        for (;count < numberOfBins; ++count)
        {
            CharString fastaFile;
            appendFileName(fastaFile, directoryPath, count);
            append(fastaFile, ext);
            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(fastaFile)))
                break;
        }
        if (count == numberOfBins)
            return ext;
    }

    // no common extensionfor all bins found.
    std::cerr << "The given directory:\n\t" << directoryPath
    << "\ndoes not contain the fasta files of all the bins!"
    << "\nAll files should have identical valid extension!\n";
    exit(1);
    return "";
}

// ----------------------------------------------------------------------------
// Function verifyFnaDir()
// ----------------------------------------------------------------------------
inline bool verifyFnaDir(CharString const directoryPath, uint32_t const numberOfBins)
{
    for (uint32_t i=0; i < numberOfBins; ++i)
    {
        CharString fastaFile;
        appendFileName(fastaFile, directoryPath, i);
        append(fastaFile, ".fna");
        if(!verifyFnaFile(fastaFile))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function appendTrailingSlash()
// ----------------------------------------------------------------------------
inline void appendTrailingSlash(CharString & directoryPath)
{
    char lastChar = directoryPath[length(directoryPath) - 1];
    if (lastChar != '/')
        append(directoryPath, "/");
}

// ----------------------------------------------------------------------------
// Function checkOutputFile()
// ----------------------------------------------------------------------------
inline bool checkOutputFile(CharString const filePath)
{
    if (std::ifstream(toCString(filePath)))
    {
        std::cerr << "File:\n" << filePath << "\nalready exists!" << std::endl;
        std::cerr << "You can either remove it or use a different file name!" << std::endl;
        return false;
    }
    std::ofstream file(toCString(filePath));
    if (!file)
    {
        std::cerr << "File:\n" << filePath << "\ncould not be created!" << std::endl;
        return false;
    }
    std::remove(toCString(filePath)); // delete file

    return true;
}

// ----------------------------------------------------------------------------
// Function getBinNoFromFile()
// ----------------------------------------------------------------------------
bool getBinNoFromFile(uint32_t & binNo, CharString const & currentFile)
{
    CharString fileName = getFilename(currentFile);
    CharString fileNameNoext = trimExtension(fileName);
    //    std::string _binNo = toCString(fileNameNoext);

    char* p;
    binNo = std::strtol(toCString(fileNameNoext), &p, 10);
    return *p == 0;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
    char comma[3] = {'\0', ' ', '\0'};
    for (const auto& e : v) {
        s << comma << e;
        comma[0] = ',';
    }
    return s;
}

#endif  // #ifndef APP_YARA_MISC_OPTIONS_DIS_H_
