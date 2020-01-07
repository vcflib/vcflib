VCF_LIB_LOCAL:=$(shell pwd)
BUILD_DIR:=$(VCF_LIB_LOCAL)/build
BIN_DIR:=$(VCF_LIB_LOCAL)/bin
LIB_DIR:=$(VCF_LIB_LOCAL)/lib
INC_DIR:=$(VCF_LIB_LOCAL)/include
CMAKE_FLAGS?=

all: 
	if [ ! -d $(BUILD_DIR) ]; then mkdir -p $(BUILD_DIR); fi
	cd $(BUILD_DIR); \
		cmake $(CMAKE_FLAGS) -DCMAKE_INSTALL_PREFIX=$(VCF_LIB_LOCAL) $(VCF_LIB_LOCAL); \
		$(MAKE) && $(MAKE) install

openmp:
	CMAKE_FLAGS=-DOPENMP=ON $(MAKE) all

profiling:
	CMAKE_FLAGS=-DPROFILING=ON $(MAKE) all

gprof:
	CMAKE_FLAGS=-DGPROF=ON $(MAKE) all

test: $(BINS)
	@prove -Itests/lib -w tests/*.t

pull:
	git pull

update: pull all

clean:
	rm -f $(BINS) $(OBJECTS)
	rm -rf $(BIN_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(INC_DIR)

.PHONY: clean all test pre
