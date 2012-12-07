cmake_minimum_required(VERSION 2.8)

add_library(vcflib
	Variant.cpp Variant.h
	split.cpp split.h convert.h join.h
)

add_subdirectory(smithwaterman)
add_subdirectory(fastahack)
add_subdirectory(multichoose)
add_subdirectory(tabixpp)

target_link_libraries(vcflib
	smithwaterman
	fastahack
	multichoose
	tabixpp
)

find_package( ZLIB REQUIRED )

add_executable(vcf2tsv
	vcf2tsv.cpp
)

target_link_libraries(vcf2tsv
	vcflib
)

add_executable(vcfaddinfo
	vcfaddinfo.cpp
)

target_link_libraries(vcfaddinfo
	vcflib
)

add_executable(vcfafpath
	vcfafpath.cpp
)

target_link_libraries(vcfafpath
	vcflib
)

add_executable(vcfallelicprimitives
	vcfallelicprimitives.cpp
)

target_link_libraries(vcfallelicprimitives
	vcflib
)

add_executable(vcfaltcount
	vcfaltcount.cpp
)

target_link_libraries(vcfaltcount
	vcflib
)


add_executable(vcfannotate
	vcfannotate.cpp
	BedReader.h
)

#add_custom_library(bedreader
#	BedReader.h
#)

target_link_libraries(vcfannotate
	bedreader
	vcflib
)

add_executable(vcfannotategenotypes
	vcfannotategenotypes.cpp
)

target_link_libraries(vcfannotategenotypes
	vcflib
)

add_executable(vcfbreakmulti
	vcfbreakmulti.cpp
)

target_link_libraries(vcfbreakmulti
	vcflib
)

add_executable(vcfclassify
	vcfclassify.cpp
)

target_link_libraries(vcfclassify
	vcflib
)

add_executable(vcfcleancomplex
	vcfcleancomplex.cpp
)

target_link_libraries(vcfcleancomplex
	vcflib
)

add_executable(vcfcommonsamples
	vcfcommonsamples.cpp
)

target_link_libraries(vcfcommonsamples
	vcflib
)

add_executable(vcfcountalleles
	vcfcountalleles.cpp
)

target_link_libraries(vcfcountalleles
	vcflib
)

add_executable(vcfdistance
	vcfdistance.cpp
)

target_link_libraries(vcfdistance
	vcflib
)

add_executable(vcfecho
	vcfecho.cpp
)

target_link_libraries(vcfecho
	vcflib
)

add_executable(vcfentropy
	vcfentropy.cpp
)

target_link_libraries(vcfentropy
	vcflib
)

add_executable(vcffilter
	vcffilter.cpp
)

target_link_libraries(vcffilter
	vcflib
)

add_executable(vcffixup
	vcffixup.cpp
)

target_link_libraries(vcffixup
	vcflib
)

add_executable(vcfflatten
	vcfflatten.cpp
)

target_link_libraries(vcfflatten
	vcflib
)

add_executable(vcfgeno2haplo
	vcfgeno2haplo.cpp
)

target_link_libraries(vcfgeno2haplo
	vcflib
)

add_executable(vcfgenotypecompare
	vcfgenotypecompare.cpp
)

target_link_libraries(vcfgenotypecompare
	vcflib
)

add_executable(vcfgenotypes
	vcfgenotypes.cpp
)

target_link_libraries(vcfgenotypes
	vcflib
)

add_executable(vcfglxgt
	vcfglxgt.cpp
)

target_link_libraries(vcfglxgt
	vcflib
)

add_executable(vcfhetcount
	vcfhetcount.cpp
)

target_link_libraries(vcfhetcount
	vcflib
)

add_executable(vcfhethomratio
	vcfhethomratio.cpp
)

target_link_libraries(vcfhethomratio
	vcflib
)

add_executable(vcfintersect
	vcfintersect.cpp
)

target_link_libraries(vcfintersect
	vcflib
)

add_executable(vcfkeepgeno
	vcfkeepgeno.cpp
)

target_link_libraries(vcfkeepgeno
	vcflib
)

add_executable(vcfkeepinfo
	vcfkeepinfo.cpp
)

target_link_libraries(vcfkeepinfo
	vcflib
)

add_executable(vcfkeepsamples
	vcfkeepsamples.cpp
)

target_link_libraries(vcfkeepsamples
	vcflib
)

add_executable(vcflength
	vcflength.cpp
)

target_link_libraries(vcflength
	vcflib
)

add_executable(vcfnumalt
	vcfnumalt.cpp
)

target_link_libraries(vcfnumalt
	vcflib
)

add_executable(vcfoverlay
	vcfoverlay.cpp
)

target_link_libraries(vcfoverlay
	vcflib
)

add_executable(vcfparsealts
	vcfparsealts.cpp
)

target_link_libraries(vcfparsealts
	vcflib
)

add_executable(vcfprimers
	vcfprimers.cpp
)

target_link_libraries(vcfprimers
	vcflib
)

add_executable(vcfrandom
	vcfrandom.cpp
)

target_link_libraries(vcfrandom
	vcflib
)

add_executable(vcfrandomsample
	vcfrandomsample.cpp
)

target_link_libraries(vcfrandomsample
	vcflib
)

add_executable(vcfremap
	vcfremap.cpp
)

target_link_libraries(vcfremap
	vcflib
)

add_executable(vcfremoveaberrantgenotypes
	vcfremoveaberrantgenotypes.cpp
)

target_link_libraries(vcfremoveaberrantgenotypes
	vcflib
)

add_executable(vcfremovesamples
	vcfremovesamples.cpp
)

target_link_libraries(vcfremovesamples
	vcflib
)

add_executable(vcfsamplediff
	vcfsamplediff.cpp
)

target_link_libraries(vcfsamplediff
	vcflib
)

add_executable(vcfsamplenames
	vcfsamplenames.cpp
)

target_link_libraries(vcfsamplenames
	vcflib
)

add_executable(vcfsamplestats
	vcfsamplestats.cpp
)

target_link_libraries(vcfsamplestats
	vcflib
)

add_executable(vcfsitesummarize
	vcfsitesummarize.cpp
)

target_link_libraries(vcfsitesummarize
	vcflib
)

add_executable(vcfstats
	vcfstats.cpp
)

target_link_libraries(vcfstats
	vcflib
)

add_executable(vcfstreamsort
	vcfstreamsort.cpp
)

target_link_libraries(vcfstreamsort
	vcflib
)

