cmake_minimum_required(VERSION 2.8)
project(vcflib)

include(ExternalProject)
include(FeatureSummary)

find_package(PkgConfig REQUIRED)

include(FindBZip2)
include(FindLibLZMA)
include(FindZLIB)

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

# ---- Options

option(BUILD_DOC "Build documentation" ON)
option(OPENMP "Enable OpenMP" OFF)
option(PROFILING "Enable profiling" OFF)
option(GPROF "Enable gprof profiling" OFF)

# ---- Dependencies

find_package(ZLIB)
set_package_properties(ZLIB PROPERTIES TYPE REQUIRED)
#find_package(LIBLZMA)
#set_package_properties(LIBLZMA PROPERTIES TYPE REQUIRED)
find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

pkg_check_modules(HTSLIB REQUIRED htslib)
pkg_check_modules(TABIXPP tabixpp) # Optional

# ---- Build switches

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++0x -D_FILE_OFFSET_BITS=64 -Isrc -O3")

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

# ---- Include files

include_directories(include)
# include_directories(tabixpp)
include_directories(fastahack)
include_directories(intervaltree)
include_directories(smithwaterman)
include_directories(multichoose)
include_directories(filevercmp)

if(HTSLIB_LOCAL)
  include_directories(tabixpp/htslib)
endif(HTSLIB_LOCAL)

file(GLOB INCLUDES
  src/*.h*
  multichoose/*.h*
  intervaltree/*.h*
  # tabixpp/*.h*
  smithwaterman/*.h*
  fastahack/*.h*
  filevercmp/*.h*)

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
    # tabixpp/tabix.cpp # do we need this?
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

# ---- Get version

find_program(GIT git)

if (GIT)

  execute_process(COMMAND git describe --abbrev=4 --dirty --always
    OUTPUT_VARIABLE GIT_VERSION
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  add_definitions(-DVERSION="${GIT_VERSION}")

else (GIT)

  file (STRINGS "VERSION" BUILD_NUMBER)
  add_definitions(-DVCFLIB_VERSION="${BUILD_NUMBER}")

endif (GIT)

# ---- Build htslib
#
# Note by default we use the distributed htslib! These are
# the old instructions:

if (HTSLIB_LOCAL)

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
# TODO: Stop vendoring in htslib and just use find_package
set(HTSLIB_EXTRA_LIBS "-lcurl" CACHE STRING "Library flags needed to link with htslib's dependencies, for chosen configuration")
set_property(TARGET htslib PROPERTY INTERFACE_LINK_LIBRARIES ${HTSLIB_EXTRA_LIBS})

endif(HTSLIB_LOCAL)

# ---- Build all

if (NOT BUILD_ONLY_LIB)
  foreach(BIN ${BINS})
    add_executable(${BIN} src/${BIN}.cpp)
    target_include_directories(${BIN} PUBLIC ${TABIXPP_LIBRARIES} ${HTSLIB_INCLUDE_DIRS})
    target_link_libraries(${BIN} vcflib tabixpp -DGIT_VERSION=0.9.0 ${TABIXPP_LIBRARIES} ${HTSLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES})
  endforeach(BIN ${BINS})
  install(TARGETS ${BINS} RUNTIME DESTINATION bin)
endif()

enable_testing()

# ---- Test

add_test(
    NAME RunTests
    COMMAND "make"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
)

function(add_pytest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 pytest/${TEST_FILE}.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

function(add_pydoctest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE pytest/${TEST_FILE}.md
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

function(add_pydoctest_fullname TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF ${TEST_FILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

# add_pytest(vcf2tsv-test)
add_pydoctest(vcf2tsv)
add_pydoctest_fullname(../README.md)

# ---- Build docs
#
# Generates man pages for the python doctests. Don't need
# to run every time so it is a separate command. For pandoc logic see
# https://www.howtogeek.com/682871/how-to-create-a-man-page-on-linux/
#
# cmake --build . --target doc ; cmake --install .

find_program(PANDOC pandoc)

if (PANDOC)
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc ALL
        COMMAND ${PANDOC} ./test/pytest/vcf2tsv.md -s -t man -o ./man/vcf2tsv.1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating man documentation with PANDOC"
        VERBATIM )
else (PANDOC)
  message("Pandoc needs to be installed to generate the man pages")
endif (PANDOC)

# ---- Install

install(TARGETS vcflib ARCHIVE DESTINATION lib)

install(FILES ${INCLUDES} DESTINATION include)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/man DESTINATION ${CMAKE_INSTALL_PREFIX}/man/man1)
