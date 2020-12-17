######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
CFLAGS=-Wall -std=c++0x
LINKERS=-lm -lpthread -lhts -lvcflib -lz
INCLUDE=-I../src -I../googletest/googletest/include/ -I../fastahack -I../tabixpp/ -I../tabixpp/htslib/ -I ../smithwaterman/ -I ../multichoose/ -I ../filevercmp/ -I ../googletest/googletest/make/
LIB=-L../build -L../tabixpp/htslib -L../tabixpp/htslib/ -L../googletest/googletest/make/ -lcurl -llzma -lbz2
LIBGTEST=../googletest/googletest/make/gtest_main.a
LIBVCF=../build/libvcflib.a
all: run

../googletest/googletest/make/gtest_main.a:
	cd ../googletest/googletest/make/ && make
tests/main: ../googletest/googletest/make/gtest_main.a
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIB) $(LINKERS) tests/mainTest.cpp ../tabixpp/tabix.cpp ../tabixpp/htslib/libhts.a $(LIBVCF) $(LIBGTEST) -o tests/main
run: tests/main
	./tests/main
clean:
	rm -f tests/main
