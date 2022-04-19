cmake_minimum_required(VERSION 3.16)
project(vcflib)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
include(FeatureSummary)

find_package(PkgConfig REQUIRED)

include(FindBZip2)
include(FindLibLZMA)
include(FindZLIB)
include(FindCURL) # for htslib

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

# ---- Options

option(BUILD_DOC "Build documentation" ON)
option(OPENMP "Enable OpenMP" ON)
option(PROFILING "Enable profiling" OFF)
option(GPROF "Enable gprof profiling" OFF)
option(ASAN "Use address sanitiser" OFF)

include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

# ---- Dependencies

find_package(ZLIB)
set_package_properties(ZLIB PROPERTIES TYPE REQUIRED)
#find_package(LIBLZMA)
#set_package_properties(LIBLZMA PROPERTIES TYPE REQUIRED)
find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

pkg_check_modules(HTSLIB htslib)   # Optionally builds from contrib/
pkg_check_modules(TABIXPP tabixpp) # Optionally builds from contrib/

# ---- Build switches
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  set(EXTRA_FLAGS "-O3 -march=native -D_FILE_OFFSET_BITS=64")
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif()

if ((${CMAKE_BUILD_TYPE} MATCHES Release) OR (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo))
  set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${EXTRA_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}")
endif ()

if (ASAN)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
endif(ASAN)

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
include_directories(fastahack)
include_directories(intervaltree)
include_directories(smithwaterman)
include_directories(multichoose)
include_directories(filevercmp)

if(NOT HTSLIB_FOUND)
  set(HTSLIB_LOCAL contrib/tabixpp/htslib)
endif()

if (NOT TABIX_FOUND)
  include_directories(contrib/tabixpp)
  include_directories(contrib/htslib)
  set(tabixpp_SOURCE
    contrib/tabixpp/tabix.cpp
  )
endif()

set(WFA_LOCAL contrib/WFA2-lib)

file(GLOB INCLUDES
  src/*.h*
  multichoose/*.h*
  intervaltree/*.h*
  smithwaterman/*.h*
  fastahack/*.h*
  filevercmp/*.h*)

set(vcflib_SOURCE
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
    fastahack/Fasta.cpp
    smithwaterman/SmithWatermanGotoh.cpp
    smithwaterman/Repeats.cpp
    smithwaterman/IndelAllele.cpp
    smithwaterman/disorder.cpp
    smithwaterman/LeftAlign.cpp
    fsom/fsom.c
    filevercmp/filevercmp.c
)

if (NOT tabixpp_FOUND)
    list(APPEND vcflib_SOURCE ${tabixpp_SOURCE})
endif ()

add_library(vcflib STATIC
    ${vcflib_SOURCE}
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
    vcfinfosummarize
    )

set(SCRIPTS
    bed2region
    bgziptabix
    # plotBfst.R
    # plotHaplotypes.R
    # plotHapLrt.R
    # plotPfst.R
    # plot_roc.r
    # plotSmoothed.R
    # plotWCfst.R
    # plotXPEHH.R
    vcf2bed.py
    vcf2sqlite.py
    vcfbiallelic
    vcfclearid
    vcfclearinfo
    vcfcomplex
    vcffirstheader
    vcfgtcompare.sh
    vcfindelproximity
    vcfindels
    vcfjoincalls
    vcfmultiallelic
    vcfmultiway
    vcfmultiwayscripts
    vcfnobiallelicsnps
    vcfnoindels
    vcfnosnps
    vcfnulldotslashdot
    vcfplotaltdiscrepancy.r
    vcfplotaltdiscrepancy.sh
    vcfplotsitediscrepancy.r
    vcfplottstv.sh
    vcfprintaltdiscrepancy.r
    vcfprintaltdiscrepancy.sh
    vcfqualfilter
    vcfregionreduce
    vcfregionreduce_and_cut
    vcfregionreduce_pipe
    vcfregionreduce_uncompressed
    vcfremovenonATGC
    vcfsnps
    vcfsort
    vcf_strip_extra_headers
    vcfvarstats
    )

# ---- Get version

file (STRINGS "VERSION" BUILD_NUMBER)
add_definitions(-DVCFLIB_VERSION="${BUILD_NUMBER}")
add_definitions(-DVERSION="${BUILD_NUMBER}")

# ---- Build htslib
#
# Note by default we use the distributed htslib! These are
# the old instructions:

if (HTSLIB_LOCAL)

  include_directories(${HTSLIB_LOCAL})

  ExternalProject_Add(htslib-EXT
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/${HTSLIB_LOCAL}"
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_COMMAND $(MAKE) lib-static
  )
  ExternalProject_Get_property(htslib-EXT SOURCE_DIR)
  ExternalProject_Get_property(htslib-EXT INSTALL_DIR)
  set(htslib_INCLUDE "${SOURCE_DIR}")
  set(htslib_static "${SOURCE_DIR}/libhts.a")
  set(htslib_lib "${htslib_static}")
  set(HTSLIB_LIBRARIES ${htslib_lib})

  add_custom_target(htslib
    DEPENDS htslib-EXT
    VERBATIM)

  set_property(TARGET htslib PROPERTY
    IMPORTED_LOCATION ${hstlib_lib})
  add_dependencies(htslib htslib)

# If the user wants to configure our HTSlib to build with libddeflate, we need
# to make sure to link against libdeflate as a transitive dependency. To do
# that, pass -DHTSLIB_EXTRA_LIBS="-ldeflate" when configuring the project with
# cmake.
# TODO: Stop vendoring in htslib and just use find_package

  set(HTSLIB_EXTRA_LIBS "-lcurl" CACHE STRING "Library flags needed to link with htslib's dependencies, for chosen configuration")
  set_property(TARGET htslib PROPERTY INTERFACE_LINK_LIBRARIES ${HTSLIB_EXTRA_LIBS})

endif(HTSLIB_LOCAL)

if (ASAN)
  set(wfa_MAKE_ARGS BUILD_TOOLS=0 BUILD_EXAMPLES=0 setup asan lib_wfa)
else()
  set(wfa_MAKE_ARGS BUILD_TOOLS=0 BUILD_EXAMPLES=0 setup lib_wfa)
endif()

if ((${CMAKE_BUILD_TYPE} MATCHES Release) OR (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo))
  set(wfa_MAKE_ARGS ${wfa_MAKE_ARGS} release)
endif()

ExternalProject_Add(wfa-EXT
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/${WFA_LOCAL}"
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE ON
    # BUILD_COMMAND $(MAKE) BUILD_TOOLS=0 BUILD_EXAMPLES=0 setup asan lib_wfa
    BUILD_COMMAND $(MAKE) ${wfa_MAKE_ARGS}
)
ExternalProject_Get_property(wfa-EXT SOURCE_DIR)
set(WFA_INCLUDE_DIRS ${SOURCE_DIR})
set(WFA_LINK_LIBRARIES ${SOURCE_DIR}/lib/libwfacpp.a ${SOURCE_DIR}/lib/libwfa.a)
add_custom_target(wfa
    DEPENDS wfa-EXT
    VERBATIM)

include_directories(${WFA_INCLUDE_DIRS})

set(vcflib_DEPS
  CURL::libcurl
  wfa
)

set(vcflib_LIBS
  ${HTSLIB_LIBRARIES}
  ${CURL_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${LIBLZMA_LIBRARIES}
  ${BZIP2_LIBRARIES}
  ${TABIXPP_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${WFA_LINK_LIBRARIES}
)

target_link_libraries(vcflib PUBLIC ${vcflib_LIBS})
add_dependencies(vcflib ${vcflib_DEPS})

# ---- Build all

if (NOT BUILD_ONLY_LIB)
  foreach(BIN ${BINS})
    add_executable(${BIN} src/${BIN}.cpp)
    # target_include_directories(${BIN} PUBLIC ${HTSLIB_INCLUDE_DIRS})
    add_dependencies(${BIN} vcflib)
    target_link_libraries(${BIN} PUBLIC vcflib ${vcflib_LIBS})
  endforeach(BIN ${BINS})
  # add_dependencies(${BINS} vcflib)
  install(TARGETS ${BINS} RUNTIME DESTINATION bin)

  # ---- Copy scripts
  foreach(SCRIPT ${SCRIPTS})
    install(PROGRAMS ./scripts/${SCRIPT} DESTINATION bin RENAME ${SCRIPT})
  endforeach(SCRIPT ${SCRIPTS})

endif()

enable_testing()

# ---- Test

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
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF pytest/${TEST_FILE}.md
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

function(add_doctest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF ../${TEST_FILE}.md
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
add_pydoctest(vcfallelicprimitives)
add_pydoctest(vcfnulldotslashdot)
add_doctest(doc/vcfintersect)
add_doctest(doc/vcffilter)
# add_pydoctest_fullname(../README.md)

# ---- Build docs
#
# Generates man pages for the python doctests. Don't need
# to run every time so it is a separate command. For pandoc logic see
# https://www.howtogeek.com/682871/how-to-create-a-man-page-on-linux/
#
# cmake --build . --target man ; cmake --install .

find_program(PANDOC pandoc)

if (PANDOC)
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( man
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${PANDOC} ./test/pytest/vcfnulldotslashdot.md -s -t man -o ./man/vcfnulldotslashdot.1
        COMMAND ruby ./test/scripts/bin2md.rb --man test/scripts/bin2md-template.erb
        COMMAND ruby ./test/scripts/bin2md.rb --man --index
        COMMAND ruby ./test/scripts/md2man
        # regenerate to allow for URLs in markdown docs
        COMMAND ruby ./test/scripts/bin2md.rb test/scripts/bin2md-template.erb
        COMMAND ruby ./test/scripts/bin2md.rb --index
    )
else (PANDOC)
  message("Pandoc needs to be installed to generate the man pages")
endif (PANDOC)

# ---- Install

install(TARGETS vcflib ARCHIVE DESTINATION lib)

install(FILES ${INCLUDES} DESTINATION include)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/man/ DESTINATION ${CMAKE_INSTALL_PREFIX}/man/man1)
