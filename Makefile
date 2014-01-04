#OBJ_DIR = ./
HEADERS = src/Variant.h \
		  src/split.h \
		  src/join.h
SOURCES = src/Variant.cpp \
		  src/split.cpp
OBJECTS= $(SOURCES:.cpp=.o)

# TODO
#vcfstats.cpp

BIN_SOURCES = src/vcfecho.cpp \
			  src/vcfaltcount.cpp \
			  src/vcfhetcount.cpp \
			  src/vcfhethomratio.cpp \
			  src/vcffilter.cpp \
			  src/vcf2tsv.cpp \
			  src/vcfgenotypes.cpp \
			  src/vcfannotategenotypes.cpp \
			  src/vcfcommonsamples.cpp \
			  src/vcfremovesamples.cpp \
			  src/vcfkeepsamples.cpp \
			  src/vcfsamplenames.cpp \
			  src/vcfgenotypecompare.cpp \
			  src/vcffixup.cpp \
			  src/vcfclassify.cpp \
			  src/vcfsamplediff.cpp \
			  src/vcfremoveaberrantgenotypes.cpp \
			  src/vcfrandom.cpp \
			  src/vcfparsealts.cpp \
			  src/vcfstats.cpp \
			  src/vcfflatten.cpp \
			  src/vcfprimers.cpp \
			  src/vcfnumalt.cpp \
			  src/vcfcleancomplex.cpp \
			  src/vcfintersect.cpp \
			  src/vcfannotate.cpp \
			  src/vcfallelicprimitives.cpp \
			  src/vcfoverlay.cpp \
			  src/vcfaddinfo.cpp \
			  src/vcfkeepinfo.cpp \
			  src/vcfkeepgeno.cpp \
			  src/vcfafpath.cpp \
			  src/vcfcountalleles.cpp \
			  src/vcflength.cpp \
			  src/vcfdistance.cpp \
			  src/vcfrandomsample.cpp \
			  src/vcfentropy.cpp \
			  src/vcfglxgt.cpp \
			  src/vcfroc.cpp \
			  src/vcfsom.cpp \
			  src/vcfcheck.cpp \
			  src/vcfstreamsort.cpp \
			  src/vcfuniq.cpp \
			  src/vcfuniqalleles.cpp \
			  src/vcfremap.cpp \
			  src/vcf2fasta.cpp \
			  src/vcfsitesummarize.cpp \
			  src/vcfbreakmulti.cpp \
			  src/vcfcreatemulti.cpp \
			  src/vcfevenregions.cpp \
			  src/vcfcat.cpp \
			  src/vcfgenosummarize.cpp \
			  src/vcfgenosamplenames.cpp \
			  src/vcfgeno2haplo.cpp \
			  src/vcfleftalign.cpp \
			  src/vcfcombine.cpp \
			  src/vcfgeno2alleles.cpp \

#BINS = $(BIN_SOURCES:.cpp=)
BINS = $(addprefix bin/,$(notdir $(BIN_SOURCES:.cpp=)))
SHORTBINS = $(notdir $(BIN_SOURCES:.cpp=))

TABIX = tabixpp/tabix.o

FASTAHACK = fastahack/Fasta.o

SMITHWATERMAN = smithwaterman/SmithWatermanGotoh.o 

REPEATS = smithwaterman/Repeats.o

INDELALLELE = smithwaterman/IndelAllele.o

DISORDER = smithwaterman/disorder.c

LEFTALIGN = smithwaterman/LeftAlign.o

FSOM = fsom/fsom.o

FILEVERCMP = filevercmp/filevercmp.o

INCLUDES = -I. -L. -Ltabixpp/ -ltabix -lz -lm

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
#CXXFLAGS = -O2
#CXXFLAGS = -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual

SSW = src/ssw.o src/ssw_cpp.o

ssw.o: src/ssw.h
ssw_cpp.o:src/ssw_cpp.h

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(HEADERS) $(TABIX)
	$(CXX) -c -o $@ src/$(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(TABIX):
	cd tabixpp && $(MAKE)

$(SMITHWATERMAN):
	cd smithwaterman && $(MAKE)

$(DISORDER): $(SMITHWATERMAN)

$(REPEATS): $(SMITHWATERMAN)

$(LEFTALIGN): $(SMITHWATERMAN)

$(INDELALLELE): $(SMITHWATERMAN)

$(FASTAHACK):
	cd fastahack && $(MAKE)

$(FSOM):
	cd fsom && $(CXX) $(CXXFLAGS) -c fsom.c -lm

$(FILEVERCMP):
	cd filevercmp && make

$(SHORTBINS):
	$(MAKE) bin/$@

$(BINS): $(BIN_SOURCES) $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK) $(DISORDER) $(LEFTALIGN) $(INDELALLELE) $(SSW) $(FSOM) $(FILEVERCMP)
	$(CXX) $(OBJECTS) $(SMITHWATERMAN) $(REPEATS) $(DISORDER) $(LEFTALIGN) $(INDELALLELE) $(SSW) $(FASTAHACK) $(FSOM) $(FILEVERCMP) tabixpp/tabix.o tabixpp/bgzf.o src/$(notdir $@).cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

clean:
	rm -f $(BINS) $(OBJECTS)
	rm -f ssw_cpp.o ssw.o
	cd tabixpp && make clean
	cd smithwaterman && make clean
	cd fastahack && make clean

.PHONY: clean all
