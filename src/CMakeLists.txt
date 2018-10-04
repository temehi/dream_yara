# ===========================================================================
#          DREAM-Yara - a distributed version of the Yara read mapper
# ===========================================================================
# File: /src/CMakeLists.txt
#
# CMakeLists.txt file for dream_yara
# ===========================================================================

# ----------------------------------------------------------------------------
# App version
# ----------------------------------------------------------------------------

# change this after every release
set (SEQAN_APP_VERSION_MAJOR "1")
set (SEQAN_APP_VERSION_MINOR "0")
set (SEQAN_APP_VERSION_PATCH "0")

# don't change the following
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

# ----------------------------------------------------------------------------
# Dependencies (continued)
# ----------------------------------------------------------------------------

if (COMPILER_MSVC)
    message (STATUS "DREAM-Yara is currently not build on Microsoft Visual
     because it requires OpenMP 3.0.")
    return ()
endif ()

# Search SeqAn and select dependencies.
find_package(OpenMP QUIET)
find_package(ZLIB   QUIET)
find_package(BZip2  QUIET)


# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------
option (DR_YARA_NATIVE_BUILD "Architecture-specific optimizations, i.e. g++ -march=native."                      ON)
option (DR_YARA_STATIC_BUILD "Include all libraries in the binaries."                                            OFF)

if (DR_YARA_NATIVE_BUILD)
    add_definitions (-DDR_YARA_NATIVE_BUILD=1)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
    if (COMPILER_IS_INTEL)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif (COMPILER_IS_INTEL)
endif (DR_YARA_NATIVE_BUILD)

if (DR_YARA_STATIC_BUILD)
    add_definitions (-DDR_YARA_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        message (WARNING "WARNING: Builds on Mac are never fully static.")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    else (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif (APPLE)
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
endif (DR_YARA_STATIC_BUILD)

# ----------------------------------------------------------------------------
# Dependencies (continued)
# ----------------------------------------------------------------------------


# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  OpenMP 3.0 not found: building dream_yara without multi-threading.")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib not found: building dream_yara without gzipped input and BAM output.")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (STATUS "  BZip2 not found: building dream_yara without bzipped input.")
endif (NOT BZIP2_FOUND)


if (COMPILER_MSVC)
    message (STATUS "DREAM-Yara is currently not build on Microsoft Visual because it requires OpenMP 3.0.")
    return ()
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
        return ()
    endif ()
elseif (COMPILER_IS_CLANG)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
        message (FATAL_ERROR "Your Clang version is too old. Please upgrade to 3.8.0 or use GCC.")
        return()
    elseif (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        message (WARNING "Your Clang version is too old, you will not have parallelism! Upgrade to 3.8.0 or newer.")
    endif ()
elseif (COMPILER_IS_INTEL)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0.2)
        message(WARNING "Your Intel Compiler version is too old. Please upgrade to 16.0.2 or newer!")
    endif ()
endif ()

find_package(SeqAn  QUIET REQUIRED CONFIG)

if (SEQAN_VERSION_STRING VERSION_LESS "2.3.1")
    message (FATAL_ERROR "The minimum SeqAn verison required is 2.3.1!")
    return ()
endif ()

message(STATUS "These dependencies where found:")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     ZLIB       ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     BZIP2      ${BZIP2_FOUND}      ${BZIP2_VERSION_STRING}")
message(   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")
message(STATUS "The requirements were met.")


message ("\n${ColourBold}Build configuration${ColourReset}")

message (STATUS "DREAM-YARA version is: ${SEQAN_APP_VERSION}")

message(STATUS "The following options are selected for the build:")
message(   "     DR_YARA_NATIVE_BUILD      ${DR_YARA_NATIVE_BUILD}")
message(   "     DR_YARA_STATIC_BUILD      ${DR_YARA_STATIC_BUILD}")
message(STATUS "Run 'cmake -LH' to get a comment on each option.")
message(STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# ----------------------------------------------------------------------------
# Compiler specifics
# ----------------------------------------------------------------------------

if (COMPILER_IS_CLANG)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ftemplate-depth-1024")

    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla-extension")
endif (COMPILER_IS_CLANG)

if(CMAKE_COMPILER_IS_GNUCXX)
    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla")

    # parallelize parts of build even for one translation unit
    include(ProcessorCount)
    ProcessorCount(NCPU)
    if(NCPU GREATER 1)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -flto=${NCPU}")
    endif()
endif(CMAKE_COMPILER_IS_GNUCXX)

if (NOT COMPILER_IS_INTEL)
    # strip binaries to make them smaller
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif ()

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "1.0.2")

option (DR_YARA_LARGE_CONTIGS "Set to OFF to disable support for more than 32k contigs or contigs longer than 4Gbp." ON)
if (DR_YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)
  add_definitions (-DDR_YARA_LARGE_CONTIGS=1)
endif (DR_YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)

# ----------------------------------------------------------------------------
# Disable Compiler Warnings
# ----------------------------------------------------------------------------

# Disables warning about conversion from 64 bit to 32 bit integrals.
# Note that there is a runtime exception thrown if the type is not sufficient.
# Unfortunately the VC does not detect it.
if (MSVC)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4244")  # When compiling on 32 bit windows.
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4267")  # When compiling on 64 bit windows.
endif (MSVC)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------
# Enable global exception handler for all seqan apps.
set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1)

# Add include directories.
include_directories (${SDSL_INCLUDE_DIRS})
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Set the right output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(dream_yara_indexer d_indexer.cpp
                                  misc_tags.h
                                  misc_timer.h
                                  misc_options.h
                                  d_misc_options.h
                                  misc_types.h
                                  store_seqs.h
                                  bits_matches.h
                                  index_fm.h)

add_executable(dream_yara_build_filter d_build_filter.cpp
                                       misc_timer.h
                                       misc_options.h
                                       d_misc_options.h
                                       misc_types.h
                                       store_seqs.h
                                       bits_matches.h
                                       d_kdx_filter.h
                                       d_bloom_filter.h)

add_executable(dream_yara_update_filter d_update_filter.cpp
                                        misc_timer.h
                                        misc_options.h
                                        d_misc_options.h
                                        misc_types.h
                                        store_seqs.h
                                        bits_matches.h
                                        d_kdx_filter.h
                                        d_bloom_filter.h)

add_executable(dream_yara_mapper  d_mapper.cpp
                                  d_mapper.h
                                  mapper.h
                                  mapper_collector.h
                                  mapper_classifier.h
                                  mapper_ranker.h
                                  mapper_filter.h
                                  mapper_extender.h
                                  mapper_verifier.h
                                  mapper_aligner.h
                                  mapper_writer.h
                                  basic_alphabet.h
                                  file_pair.h
                                  file_prefetched.h
                                  store_seqs.h
                                  misc_tags.h
                                  misc_timer.h
                                  misc_options.h
                                  d_misc_options.h
                                  d_kdx_filter.h
                                  d_bloom_filter.h
                                  misc_types.h
                                  bits_reads.h
                                  bits_hits.h
                                  bits_matches.h
                                  bits_context.h
                                  bits_seeds.h
                                  bits_bucket.h
                                  find_extender.h
                                  find_verifier.h
                                  index_fm.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (dream_yara_indexer ${SEQAN_LIBRARIES})
target_link_libraries (dream_yara_build_filter ${SEQAN_LIBRARIES})
target_link_libraries (dream_yara_update_filter ${SEQAN_LIBRARIES})
target_link_libraries (dream_yara_mapper ${SEQAN_LIBRARIES})


# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------
# Install dream_yara in ${PREFIX}/bin directory
install (TARGETS  dream_yara_indexer
                  dream_yara_build_filter
                  dream_yara_update_filter
                  dream_yara_mapper
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/dream_yara for SeqAn release builds.
install (FILES LICENSE README.rst
         DESTINATION "share/doc/dream_yara")


# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable yara_mapper in CTD structure.
set (SEQAN_CTD_EXECUTABLES
                  ${SEQAN_CTD_EXECUTABLES}
                  dream_yara_indexer
                  dream_yara_build_filter
                  dream_yara_update_filter
                  dream_yara_mapper
    CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "dream_yara")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "DREAM-Yara - a distributed version of the Yara read mapper")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>")
set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.rst")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
set (CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
set (CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_APP_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${SEQAN_APP_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${SEQAN_APP_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CPACK_GENERATOR "ZIP;NSIS")
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CPACK_GENERATOR "ZIP;DragNDrop")
elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
    set(CPACK_GENERATOR "TBZ2")
else()
    set(CPACK_GENERATOR "TXZ")
endif ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
endif ()

# Package architecture
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
    set(CMAKE_SYSTEM_PROCESSOR "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
endif ()
if (CMAKE_CXX_FLAGS MATCHES "avx2")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_avx2")
elseif (CMAKE_CXX_FLAGS MATCHES "sse4")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_sse4")
endif()

# Include architecture in package name
if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

#replace Darwin with Mac
string(REPLACE "Darwin" "Mac" CPACK_SYSTEM_NAME ${CPACK_SYSTEM_NAME})

include (CPack)
