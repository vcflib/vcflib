vcfstats: Variant.h Split.cpp Split.h vcfstats.cpp
	g++ -O3 vcfstats.cpp Split.cpp -o vcfstats

clean:
	rm -f vcfstats *.o

.PHONY: clean
