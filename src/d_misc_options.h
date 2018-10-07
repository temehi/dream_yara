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

static const uint32_t BD_METADATA_SIZE = 256;
static const uint8_t INT_WIDTH = 0x40;


enum FilterType
{
    BLOOM, KMER_DIRECT, NONE
};

// ============================================================================
// Functions
// ============================================================================
// ----------------------------------------------------------------------------
// Function append_file_name()
// ----------------------------------------------------------------------------
inline void append_file_name(CharString & target, CharString const & source, uint32_t const i)
{
    target = source;
    append(target, std::to_string(i));
}

inline void append_file_name(CharString & target, uint32_t const i)
{
    CharString source = target;
    append_file_name(target, source, i);
}

// ----------------------------------------------------------------------------
// Function get_ext_with_leading_dot()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString const>::Type
get_ext_with_leading_dot(TString const & string)
{
    return suffix(string, firstOf(string, IsDot()));
}

// ----------------------------------------------------------------------------
// Function get_files_in_dir()
// ----------------------------------------------------------------------------
inline void get_files_in_dir(StringSet<CharString> & file_names, CharString const directory_path)
{
    DIR *dir;
    struct dirent *ent;
    struct stat st;

    dir = opendir(toCString(directory_path));
    while ((ent = readdir(dir)) != NULL)
    {
        CharString file_name = ent->d_name;
        CharString full_file_name = directory_path;
        append(full_file_name, "/");
        append(full_file_name, file_name);


        bool invalidFile = (file_name[0] == '.')
        || (stat(toCString(full_file_name), &st) == -1)
        || ((st.st_mode & S_IFDIR) != 0);

        if (!invalidFile)
            appendValue(file_names, file_name);
    }
}

// ----------------------------------------------------------------------------
// Function get_valid_files_in_dir()
// ----------------------------------------------------------------------------
inline void get_valid_files_in_dir(StringSet<CharString> & file_names,
                                   CharString const directory_path,
                                   std::vector<std::string> const & valid_ext)
{
    StringSet<CharString>  all_file_names;
    get_files_in_dir(all_file_names, directory_path);
    for (uint32_t i = 0; i < length(all_file_names); ++i)
    {
        CharString ext = get_ext_with_leading_dot(all_file_names[i]);

        auto it = std::find(valid_ext.begin(), valid_ext.end(), std::string(toCString(ext)));
        if (it != valid_ext.end())
            appendValue(file_names, all_file_names[i]);
    }

}

// ----------------------------------------------------------------------------
// Function verify_indices_dir()
// ----------------------------------------------------------------------------
inline bool verify_indices_dir(CharString const directory_path, uint32_t const number_of_bins)
{
    for (uint32_t i=0; i < number_of_bins; ++i)
    {
        CharString contigs_limit_file;
        append_file_name(contigs_limit_file, directory_path, i);
        append(contigs_limit_file, ".txt.size");

        String<uint64_t> limits;

        if (!open(limits, toCString(contigs_limit_file), OPEN_RDONLY|OPEN_QUIET))
        {
            std::cerr << "No index for bin " << i << '\n';
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function verify_fna_file()
// ----------------------------------------------------------------------------
inline bool verify_fna_file(CharString const & fasta_file)
{
    SeqFileIn seq_file_in;
    if (!open(seq_file_in, toCString(fasta_file)))
    {
        std::cerr << "Fasta file: " << fasta_file << " can not be found!\n" ;
        return false;
    }
    close(seq_file_in);
    return true;
}

// ----------------------------------------------------------------------------
// Function common_ext()
// ----------------------------------------------------------------------------
inline std::string common_ext(CharString const directory_path, uint32_t const number_of_bins)
{
    std::vector<std::string> extensions =  SeqFileIn::getFileExtensions();

    uint32_t count = 0;
    for (auto ext : extensions)
    {
        count = 0;
        for (;count < number_of_bins; ++count)
        {
            CharString fasta_file;
            append_file_name(fasta_file, directory_path, count);
            append(fasta_file, ext);
            SeqFileIn seq_file_in;
            if (!open(seq_file_in, toCString(fasta_file)))
                break;
        }
        if (count == number_of_bins)
            return ext;
    }

    // no common extensionfor all bins found.
    std::cerr << "The given directory:\n\t" << directory_path
    << "\ndoes not contain the fasta files of all the bins!"
    << "\nAll files should have identical valid extension!\n";
    exit(1);
    return "";
}

// ----------------------------------------------------------------------------
// Function verify_fna_dir()
// ----------------------------------------------------------------------------
inline bool verify_fna_dir(CharString const directory_path, uint32_t const number_of_bins)
{
    for (uint32_t i=0; i < number_of_bins; ++i)
    {
        CharString fasta_file;
        append_file_name(fasta_file, directory_path, i);
        append(fasta_file, ".fna");
        if(!verify_fna_file(fasta_file))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function append_trailing_slash()
// ----------------------------------------------------------------------------
inline void append_trailing_slash(CharString & directory_path)
{
    char lastChar = directory_path[length(directory_path) - 1];
    if (lastChar != '/')
        append(directory_path, "/");
}

// ----------------------------------------------------------------------------
// Function check_output_file()
// ----------------------------------------------------------------------------
inline bool check_output_file(CharString const file_path)
{
    if (std::ifstream(toCString(file_path)))
    {
        std::cerr << "File:\n" << file_path << "\nalready exists!" << std::endl;
        std::cerr << "You can either remove it or use a different file name!" << std::endl;
        return false;
    }
    std::ofstream file(toCString(file_path));
    if (!file)
    {
        std::cerr << "File:\n" << file_path << "\ncould not be created!" << std::endl;
        return false;
    }
    std::remove(toCString(file_path)); // delete file

    return true;
}
// ----------------------------------------------------------------------------
// Function get_file_name()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString const>::Type
get_file_name(TString const & string)
{
    return suffix(string, lastOf(string, IsPathDelimited()));
}

// ----------------------------------------------------------------------------
// Function get_bin_number_from_file()
// ----------------------------------------------------------------------------
bool get_bin_number_from_file(uint32_t & bin_number, CharString const & current_file)
{
    CharString file_name = get_file_name(current_file);
    CharString file_name_no_ext = trimExtension(file_name);
    //    std::string _bin_number = toCString(file_name_no_ext);

    char* p;
    bin_number = std::strtol(toCString(file_name_no_ext), &p, 10);
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
