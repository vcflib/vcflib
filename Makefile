vcfecho: Variant.h Split.cpp Split.h vcfecho.cpp
	g++ -g -O3 vcfecho.cpp Split.cpp -o vcfecho

clean:
	rm -f vcfecho *.o

.PHONY: clean
