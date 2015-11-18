#OBJ_DIR = ./
HEADERS = src/Variant.h \
		  src/split.h \
		  src/join.h
SOURCES = src/Variant.cpp \
		  src/split.cpp
OBJECTS= $(SOURCES:.cpp=.o)

VCF_LIB_LOCAL:=$(shell pwd)
BIN_DIR:=bin
LIB_DIR:=lib
SRC_DIR=src
INC_DIR:=include
OBJ_DIR:=obj

MAKE ?=		make

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
		src/vcfindex.cpp \
		src/vcf2dag.cpp \
		src/vcfsample2info.cpp \
		src/vcfqual2info.cpp \
		src/vcfinfo2qual.cpp \
		src/vcfglbound.cpp \
		src/vcfinfosummarize.cpp

# when we can figure out how to build on mac
# src/vcfsom.cpp

#BINS = $(BIN_SOURCES:.cpp=)
BINS = $(addprefix bin/,$(notdir $(BIN_SOURCES:.cpp=)))
SHORTBINS = $(notdir $(BIN_SOURCES:.cpp=))

# Use ?= to allow overriding submodule install locations from the env
# or command-line
SW_PATH ?=	smithwaterman
TABIX_PATH ?=	tabixpp
HTS_PATH ?=	$(TABIX_PATH)/htslib
FH_PATH ?=	fastahack

# FIXME: Replace each of these with a library, like SMITHWATERMAN
# and allow overriding the prefix
SMITHWATERMAN ?=	$(SW_PATH)/libsw.a
TABIX ?=		$(TABIX_PATH)/libtabix.a
FASTAHACK ?=		$(FH_PATH)/libfastahack.a
FSOM =			fsom/fsom.o
FILEVERCMP =		filevercmp/filevercmp.o

INCLUDES =	-I$(HTS_PATH) -I$(INC_DIR)
LDFLAGS =	-L$(LIB_DIR) -lvcflib \
		-L$(SW_PATH) -lsw \
		-L$(TABIX_PATH) -ltabix \
		-L$(HTS_PATH) -lhts \
		-L$(FH_PATH) -lfastahack \
		-lpthread -lz -lm

all: $(OBJECTS) $(BINS)

CXX ?=		g++
CXXFLAGS ?=	-O3
#CXXFLAGS ?=	-O2
CXXFLAGS +=	-D_FILE_OFFSET_BITS=64
#CXXFLAGS +=	-pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual

SSW = src/ssw.o src/ssw_cpp.o

ssw.o: src/ssw.h
ssw_cpp.o:src/ssw_cpp.h

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

SUBMOD_OBJS ?= $(TABIX) multichoose $(SMITHWATERMAN) $(FILEVERCMP)

$(OBJECTS): $(SOURCES) $(HEADERS) pre $(SUBMOD_OBJS)
	$(CXX) -c -o $@ src/$(*F).cpp $(INCLUDES) $(CXXFLAGS) $(LDFLAGS) && cp src/*.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

multichoose: pre
	cd multichoose && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

intervaltree: pre
	cd intervaltree && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(TABIX): pre
	cd tabixpp && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(SMITHWATERMAN): pre
	cd smithwaterman && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && cp *.o $(VCF_LIB_LOCAL)/$(OBJ_DIR)/

$(FASTAHACK): pre
	cd fastahack && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && cp Fasta.o $(VCF_LIB_LOCAL)/$(OBJ_DIR)/

#$(FSOM):
#	cd fsom && $(CXX) $(CXXFLAGS) -c fsom.c -lm

$(FILEVERCMP): pre
	cd filevercmp && $(MAKE) && cp *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && cp *.o $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(SHORTBINS): pre
	$(MAKE) bin/$@

SUBMOD_BINS ?= $(SMITHWATERMAN) $(FASTAHACK) $(FILEVERCMP) intervaltree

$(BINS): $(BIN_SOURCES) libvcflib.a $(OBJECTS) $(SSW) pre $(SUBMOD_BINS)
	$(CXX) src/$(notdir $@).cpp -o $@ $(INCLUDES) $(CXXFLAGS) $(LDFLAGS)

libvcflib.a: $(OBJECTS) $(SMITHWATERMAN) $(SSW) $(FILEVERCMP) pre
	ar rs libvcflib.a $(OBJECTS) $(SSW) $(FILEVERCMP)
	cp libvcflib.a $(LIB_DIR)


test: $(BINS)
	@prove -Itests/lib -w tests/*.t

pre:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi

clean:
	rm -f $(BINS) $(OBJECTS)
	rm -f ssw_cpp.o ssw.o
	rm -f libvcflib.a
	rm -rf $(BIN_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)
	rm -rf $(OBJ_DIR)
	cd tabixpp && $(MAKE) clean
	cd smithwaterman && $(MAKE) clean
	cd fastahack && $(MAKE) clean

.PHONY: clean all test pre
