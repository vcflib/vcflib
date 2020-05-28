cmake_minimum_required(VERSION 2.8)
project(vcflib)

option(OPENMP "Enable OpenMP" OFF)
option(PROFILING "Enable profiling" OFF)
option(GPROF "Enable gprof profiling" OFF)

include(ExternalProject)
include(FeatureSummary)
#include(FindLibLZMA)
#include(ExternalProject)
include(FindZLIB)
include(FindBZip2)
include(FindLibLZMA)

find_package(ZLIB)
set_package_properties(ZLIB PROPERTIES TYPE REQUIRED)
#find_package(LIBLZMA)
#set_package_properties(LIBLZMA PROPERTIES TYPE REQUIRED)
find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -D_FILE_OFFSET_BITS=64")

if(OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP)

if(PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(PROFILING)

if(GPROF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(GPROF)

include_directories(include)
include_directories(tabixpp)
include_directories(tabixpp/htslib)
include_directories(fastahack)
include_directories(intervaltree)
include_directories(smithwaterman)
include_directories(multichoose)
include_directories(filevercmp)

add_library(vcflib STATIC 
    src/Variant.h 
    src/split.h 
    src/pdflib.hpp 
    src/var.hpp 
    src/cdflib.hpp 
    src/rnglib.hpp 
    src/join.h
    src/Variant.cpp 
    src/rnglib.cpp 
    src/var.cpp 
    src/pdflib.cpp 
    src/cdflib.cpp 
    src/split.cpp
    src/ssw.hpp
    src/ssw_cpp.hpp
    tabixpp/tabix.cpp
    fastahack/Fasta.cpp
    smithwaterman/SmithWatermanGotoh.cpp
    smithwaterman/Repeats.cpp
    smithwaterman/IndelAllele.cpp
    smithwaterman/disorder.cpp
    smithwaterman/LeftAlign.cpp
    fsom/fsom.c
    filevercmp/filevercmp.c
    )

set(BINS 
    vcfecho
    dumpContigsFromHeader
    bFst
    pVst
    hapLrt
    popStats
    wcFst
    iHS
    segmentFst
    segmentIhs
    genotypeSummary
    sequenceDiversity
    pFst
    smoother
    vcfld
    plotHaps
    abba-baba
    permuteGPAT++
    permuteSmooth
    normalize-iHS
    meltEHH
    vcfaltcount
    vcfhetcount
    vcfhethomratio
    vcffilter
    vcf2tsv
    vcfgenotypes
    vcfannotategenotypes
    vcfcommonsamples
    vcfremovesamples
    vcfkeepsamples
    vcfsamplenames
    vcfgenotypecompare
    vcffixup
    vcfclassify
    vcfsamplediff
    vcfremoveaberrantgenotypes
    vcfrandom
    vcfparsealts
    vcfstats
    vcfflatten
    vcfprimers
    vcfnumalt
    vcfcleancomplex
    vcfintersect
    vcfannotate
    vcfallelicprimitives
    vcfoverlay
    vcfaddinfo
    vcfkeepinfo
    vcfkeepgeno
    vcfafpath
    vcfcountalleles
    vcflength
    vcfdistance
    vcfrandomsample
    vcfentropy
    vcfglxgt
    vcfroc
    vcfcheck
    vcfstreamsort
    vcfuniq
    vcfuniqalleles
    vcfremap
    vcf2fasta
    vcfsitesummarize
    vcfbreakmulti
    vcfcreatemulti
    vcfevenregions
    vcfcat
    vcfgenosummarize
    vcfgenosamplenames
    vcfgeno2haplo
    vcfleftalign
    vcfcombine
    vcfgeno2alleles
    vcfindex
    vcf2dag
    vcfsample2info
    vcfqual2info
    vcfinfo2qual
    vcfglbound
    vcfinfosummarize)

execute_process(COMMAND git describe --abbrev=4 --dirty --always 
  OUTPUT_VARIABLE GIT_VERSION 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_STRIP_TRAILING_WHITESPACE)

ExternalProject_Add(htslib-build
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/tabixpp/htslib
  BINARY_DIR ${CMAKE_SOURCE_DIR}/tabixpp/htslib
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) lib-static
  INSTALL_COMMAND ""
  )

add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib PROPERTY 
  IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/tabixpp/htslib/libhts.a)
add_dependencies(htslib htslib-build)
# If the user wants to configure our HTSlib to build with libddeflate, we need
# to make sure to link against libdeflate as a transitive dependency. To do
# that, pass -DHTSLIB_EXTRA_LIBS="-ldeflate" when configuring the project with
# cmake.
set_property(TARGET htslib PROPERTY INTERFACE_LINK_LIBRARIES ${HTSLIB_EXTRA_LIBS})

add_definitions(-DVERSION="${GIT_VERSION}")

if (NOT BUILD_ONLY_LIB)
  foreach(BIN ${BINS})
    add_executable(${BIN} src/${BIN}.cpp)
    target_link_libraries(${BIN} vcflib htslib ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES})
  endforeach(BIN ${BINS})
  install(TARGETS ${BINS} RUNTIME DESTINATION bin)
endif()

install(TARGETS vcflib ARCHIVE DESTINATION lib)

file(GLOB INCLUDES
  src/*.h*
  multichoose/*.h*
  intervaltree/*.h*
  tabixpp/*.h*
  smithwaterman/*.h*
  fastahack/*.h*
  filevercmp/*.h*)

install(FILES ${INCLUDES} DESTINATION include)

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES 
  WHAT REQUIRED_PACKAGES_NOT_FOUND)
