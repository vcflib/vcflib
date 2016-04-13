VCF_LIB_LOCAL:=$(shell pwd)
BUILD_DIR:=$(VCF_LIB_LOCAL)/build
BIN_DIR:=$(VCF_LIB_LOCAL)/bin
LIB_DIR:=$(VCF_LIB_LOCAL)/lib
SRC_DIR=$(VCF_LIB_LOCAL)/src
INC_DIR:=$(VCF_LIB_LOCAL)/include
OBJ_DIR:=$(VCF_LIB_LOCAL)/obj

# TODO

all: 
	if [ ! -d $(BUILD_DIR) ]; then mkdir -p $(BUILD_DIR); fi
	cd $(BUILD_DIR); cmake -DCMAKE_INSTALL_PREFIX=$(VCF_LIB_LOCAL) $(VCF_LIB_LOCAL); $(MAKE) && $(MAKE) install

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

test: $(BINS)
	@prove -Itests/lib -w tests/*.t

pull:
	git pull

update: pull all

clean:
	rm -f $(BINS) $(OBJECTS)
	rm -f ssw_cpp.o ssw.o
	rm -f libvcflib.a
	rm -rf $(BIN_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)
	rm -rf $(OBJ_DIR)
	cd tabixpp && make clean
	cd smithwaterman && make clean
	cd fastahack && make clean

.PHONY: clean all test pre
