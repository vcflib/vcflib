######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
CFLAGS=-Wall -std=c++0x
LINKERS=-lm -lpthread -lhts -lvcflib -lz
INCLUDE=-I../src -I../googletest/googletest/include/ -I../fastahack -I ../smithwaterman/ -I ../multichoose/ -I ../filevercmp/ -I ../googletest/googletest/make/
LIB=-L../build -L../googletest/googletest/make/ -lhts -lcurl -llzma -lbz2
LIBGTEST=../googletest/googletest/make/gtest_main.a
LIBVCF=../build/libvcflib.a
all: run

../googletest/googletest/make/gtest_main.a:
	cd ../googletest/googletest/make/ && make

tests/main: ../googletest/googletest/make/gtest_main.a
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIB) $(LINKERS) tests/mainTest.cpp ../tabixpp/tabix.cpp $(LIBVCF) $(LIBGTEST) -o tests/main

run: tests/main
	./tests/main
	rm -v ../googletest/googletest/make/gtest_main.a
	rm -v ../googletest/googletest/make/sample1_unittest

clean:
	rm -f tests/main
