#OBJ_DIR = ./
HEADERS = src/Variant.h \
		  src/split.h \
		  src/pdflib.hpp \
		  src/var.hpp \
                  src/cdflib.hpp \
		  src/rnglib.hpp \
		  src/join.h
SOURCES = src/Variant.cpp \
		  src/rnglib.cpp \
		  src/var.cpp \
		  src/pdflib.cpp \
		  src/cdflib.cpp \
		  src/split.cpp
OBJECTS= $(SOURCES:.cpp=.o)

VCF_LIB_LOCAL:=$(shell pwd)
BIN_DIR:=bin
LIB_DIR:=lib
SRC_DIR=src
INC_DIR:=include
OBJ_DIR:=obj

RM = rm -f
CP = cp -f

# TODO
#vcfstats.cpp

BIN_SOURCES = src/vcfecho.cpp \
			  src/vcfnormalizesvs.cpp \
			  src/dumpContigsFromHeader.cpp \
			  src/bFst.cpp \
			  src/pVst.cpp \
			  src/hapLrt.cpp \
			  src/popStats.cpp \
			  src/wcFst.cpp \
			  src/iHS.cpp \
			  src/segmentFst.cpp \
			  src/segmentIhs.cpp \
			  src/genotypeSummary.cpp \
			  src/sequenceDiversity.cpp \
			  src/pFst.cpp \
			  src/smoother.cpp \
			  src/vcfld.cpp \
			  src/plotHaps.cpp \
			  src/abba-baba.cpp \
			  src/permuteGPAT++.cpp \
			  src/permuteSmooth.cpp \
			  src/normalize-iHS.cpp \
			  src/meltEHH.cpp \
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
			  src/vcfunphase.cpp \
			  src/vcfnull2ref.cpp \
			  src/vcfinfosummarize.cpp

# when we can figure out how to build on mac
# src/vcfsom.cpp

#BINS = $(BIN_SOURCES:.cpp=)
BINS = $(addprefix $(BIN_DIR)/,$(notdir $(BIN_SOURCES:.cpp=)))
SHORTBINS = $(notdir $(BIN_SOURCES:.cpp=))

TABIX = tabixpp/tabix.o
FASTAHACK = fastahack/Fasta.o
SMITHWATERMAN = smithwaterman/SmithWatermanGotoh.o
REPEATS = smithwaterman/Repeats.o
INDELALLELE = smithwaterman/IndelAllele.o
DISORDER = smithwaterman/disorder.o
LEFTALIGN = smithwaterman/LeftAlign.o
FSOM = fsom/fsom.o
FILEVERCMP = filevercmp/filevercmp.o

# Work out how to find htslib
# Use the one we ship in tabixpp unless told otherwise by the environment
HTS_LIB ?= $(VCF_LIB_LOCAL)/tabixpp/htslib/libhts.a
HTS_INCLUDES ?= -I$(VCF_LIB_LOCAL)/tabixpp/htslib
HTS_LDFLAGS ?= -L$(VCF_LIB_LOCAL)/tabixpp/htslib -lhts -lbz2 -lm -lz -llzma -pthread


INCLUDES = $(HTS_INCLUDES) -I$(INC_DIR) 
LDFLAGS = -L$(LIB_DIR) -lvcflib $(HTS_LDFLAGS) -lpthread -lz -lm -llzma -lbz2



all: $(OBJECTS) $(BINS) scriptToBin

scriptToBin: $(BINS)
	$(CP) scripts/* $(BIN_DIR)

GIT_VERSION += $(shell git describe --abbrev=4 --dirty --always)

CXXFLAGS = -Ofast -D_FILE_OFFSET_BITS=64 -std=c++0x
#CXXFLAGS = -O2
#CXXFLAGS = -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual

SSW = src/ssw.o src/ssw_cpp.o

ssw.o: src/ssw.hpp
ssw_cpp.o:src/ssw_cpp.hpp

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(HEADERS) $(TABIX) multichoose pre $(SMITHWATERMAN) $(FILEVERCMP) $(FASTAHACK)
	$(CXX) -c -o $@ src/$(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) && $(CP) src/*.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

multichoose: pre
	cd multichoose && $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

intervaltree: pre
	cd intervaltree && $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(TABIX): pre
	cd tabixpp && INCLUDES="$(HTS_INCLUDES)" LIBPATH="-L. $(HTS_LDFLAGS)" HTSLIB="$(HTS_LIB)" $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(SMITHWATERMAN): pre
	cd smithwaterman && $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && $(CP) *.o $(VCF_LIB_LOCAL)/$(OBJ_DIR)/

$(DISORDER): $(SMITHWATERMAN)

$(REPEATS): $(SMITHWATERMAN)

$(LEFTALIGN): $(SMITHWATERMAN)

$(INDELALLELE): $(SMITHWATERMAN)

$(FASTAHACK): pre
	cd fastahack && $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && $(CP) Fasta.o $(VCF_LIB_LOCAL)/$(OBJ_DIR)/

#$(FSOM):
#	cd fsom && $(CXX) $(CXXFLAGS) -c fsom.c -lm

$(FILEVERCMP): pre
	cd filevercmp && $(MAKE) && $(CP) *.h* $(VCF_LIB_LOCAL)/$(INC_DIR)/ && $(CP) *.o $(VCF_LIB_LOCAL)/$(INC_DIR)/

$(SHORTBINS): pre
	$(MAKE) $(BIN_DIR)/$@

$(BINS): $(BIN_SOURCES) libvcflib.a $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK) $(DISORDER) $(LEFTALIGN) $(INDELALLELE) $(SSW) $(FILEVERCMP) pre intervaltree
	$(CXX) src/$(notdir $@).cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -DVERSION=\"$(GIT_VERSION)\"

libvcflib.a: $(OBJECTS) $(SMITHWATERMAN) $(REPEATS) $(FASTAHACK) $(DISORDER) $(LEFTALIGN) $(INDELALLELE) $(SSW) $(FILEVERCMP) $(TABIX) pre
	ar rs libvcflib.a $(OBJECTS) smithwaterman/sw.o $(FASTAHACK) $(SSW) $(FILEVERCMP) $(TABIX)
	$(CP) libvcflib.a $(LIB_DIR)


test: $(BINS)
	@prove -Itests/lib -w tests/*.t

pre:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi


pull:
	git pull

update: pull all

clean:
	$(RM) $(BINS) $(OBJECTS)
	$(RM) ssw_cpp.o ssw.o
	$(RM) libvcflib.a
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(OBJ_DIR)
	$(MAKE) clean -C tabixpp
	$(MAKE) clean -C smithwaterman
	$(MAKE) clean -C fastahack
	$(MAKE) clean -C multichoose
	$(MAKE) clean -C fsom
	$(MAKE) clean -C libVCFH
	$(MAKE) clean -C test
	$(MAKE) clean -C filevercmp
	$(MAKE) clean -C intervaltree
	
.PHONY: clean all test pre
