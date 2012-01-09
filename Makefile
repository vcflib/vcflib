#OBJ_DIR = ./
HEADERS = Variant.h \
		  split.h \
		  join.h
SOURCES = Variant.cpp \
		  split.cpp
OBJECTS= $(SOURCES:.cpp=.o)

# TODO
#vcfstats.cpp

BIN_SOURCES = vcfecho.cpp \
			  vcfaltcount.cpp \
			  vcfhetcount.cpp \
			  vcffilter.cpp \
			  vcf2tsv.cpp \
			  vcfgenotypes.cpp \
			  vcfannotategenotypes.cpp \
			  vcfcommonsamples.cpp \
			  vcfremovesamples.cpp \
			  vcfkeepsamples.cpp \
			  vcfsamplenames.cpp \
			  vcfgenotypecompare.cpp \
			  vcffixup.cpp \
			  vcfclassify.cpp \
			  vcfsamplediff.cpp \
			  vcfremoveaberrantgenotypes.cpp \
			  vcfrandom.cpp \
			  vcfparsealts.cpp \
			  vcfstats.cpp \
			  vcfflatten.cpp \
			  vcfprimers.cpp \
			  vcfnumalt.cpp \
			  vcfcleancomplex.cpp \
			  vcfintersect.cpp \
			  vcfannotate.cpp \
			  vcfallelicprimitives.cpp \
			  vcfoverlay.cpp \
			  vcfaddinfo.cpp \
			  vcfkeepinfo.cpp \
			  vcfkeepgeno.cpp \
			  vcfafpath.cpp \
			  vcfcountalleles.cpp \
			  vcfdistance.cpp \
			  vcfrandomsample.cpp

BINS = $(BIN_SOURCES:.cpp=)

TABIX = tabixpp/tabix.o

FASTAHACK = fastahack/Fasta.o

SMITHWATERMAN = smithwaterman/SmithWatermanGotoh.o

INCLUDES = -lm -lz -L. -Ltabixpp/ -ltabix

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(HEADERS) $(TABIX)
	$(CXX) -c -o $@ $(*F).cpp  $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(TABIX):
	cd tabixpp && $(MAKE)

$(SMITHWATERMAN):
	cd smithwaterman && $(MAKE)

$(FASTAHACK):
	cd fastahack && $(MAKE)

$(BINS): $(BIN_SOURCES) $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK)
	$(CXX) $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK) tabixpp/tabix.o tabixpp/bgzf.o $@.cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd tabixpp && make clean
	cd smithwaterman && make clean
	cd fastahack && make clean

.PHONY: clean all
