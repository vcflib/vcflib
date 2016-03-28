######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
CFLAGS=-Wall -std=c++0x 
LINKERS=-lm -lpthread -lhts -lvcflib -lz 
INCLUDE=-I../src -I../googletest/googletest/include/ -I../tabixpp/ -I../tabixpp/htslib/ -I ../smithwaterman/ -I ../multichoose/ -I ../filevercmp/ -I ../googletest/googletest/make/
LIB=-L../ -L../tabixpp/htslib -L../tabixpp/htslib/ -L../googletest/googletest/make/
LIBGTEST=../googletest/googletest/make/gtest_main.a

../googletest/googletest/make/gtest_main.a:
	cd ../googletest/googletest/make/ && make
test/main:
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIB) $(LINKERS) tests/mainTest.cpp ../tabixpp/tabix.cpp ../tabixpp/htslib/libhts.a $(LIBGTEST) -o tests/main 