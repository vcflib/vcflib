#!/bin/bash

if [ $# != 3 ];
then
    echo "usage: $0 [annotation] [fileA] [fileB]"
    echo "annotates records in the first file with genotypes and sites from the second"
    exit
fi

annotation=$1
fileA=$2
fileB=$3

vcfcommonsamples $fileA $fileB \
 | vcfannotategenotypes $annotation - $fileB \
 | vcfgenotypecompare $annotation -
