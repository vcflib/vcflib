OBJ_DIR = ./
SOURCES = Variant.cpp \
		  Split.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

CXX = g++
CXXFLAGS = -O3

$(BUILT_OBJECTS): $(SOURCES)
    $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

BINS = vcfecho vcfaltcount vcfhetcount
BIN_SOURCES = vcfecho.cpp vcfaltcount.cpp vcfhetcount.cpp

$(BINS): $(BIN_SOURCES) $(BUILT_OBJECTS)
    $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

.PHONY: all

#vcfecho: Variant.o Split.o Split.h vcfecho.cpp
#	g++ -O3 vcfecho.cpp Split.o Variant.o -o vcfecho

#vcfaltcount: Variant.o Split.o Split.h vcfaltcount.cpp
#	g++ -O3 vcfaltcount.cpp Split.o Variant.o -o vcfaltcount

#vcfhetcount: Variant.o Split.o Split.h vcfhetcount.cpp
#	g++ -O3 vcfhetcount.cpp Split.o Variant.o -o vcfhetcount

#Variant.o: Variant.h Variant.cpp
#	g++ -O3 -c Variant.cpp

#Split.o: Split.h Split.cpp
#	g++ -O3 -c Split.cpp

clean:
	rm -f $(BINS) $(BUILT_OBJECTS)

.PHONY: clean all
