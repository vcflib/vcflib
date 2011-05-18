#OBJ_DIR = ./
HEADERS = Variant.h \
		  split.h \
		  join.h
SOURCES = Variant.cpp \
		  split.cpp
OBJECTS= $(SOURCES:.cpp=.o)

BIN_SOURCES = vcfecho.cpp \
			  vcfaltcount.cpp \
			  vcfhetcount.cpp \
			  vcfstats.cpp \
			  vcffilter.cpp \
			  vcf2tsv.cpp \
			  vcfgenotypes.cpp \
			  vcfannotategenotypes.cpp \
			  vcfcommonsamples.cpp \
			  vcfgenotypecompare.cpp \
			  vcffixup.cpp \
			  vcfclassify.cpp \
			  vcfsamplediff.cpp \
			  vcfremoveaberrantgenotypes.cpp

BINS = $(BIN_SOURCES:.cpp=)

TABIX = tabixpp/tabix.o

INCLUDES = -lm -lz -L. -Ltabixpp/ -ltabix

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(HEADERS) $(TABIX)
	$(CXX) -c -o $@ $(*F).cpp  $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(TABIX):
	cd tabixpp && make

$(BINS): $(BIN_SOURCES) $(OBJECTS)
	$(CXX) $(OBJECTS) tabixpp/tabix.o tabixpp/bgzf.o $@.cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd tabixpp && make clean

.PHONY: clean all
