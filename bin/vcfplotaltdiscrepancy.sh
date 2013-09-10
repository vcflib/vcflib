#!/bin/bash

filename=$1
tag=$2

vcf2tsv \
    | tsvsplit \
        QUAL \
        AC \
        $tag.has_variant \
    | tf2binary \
    | vcfplotsitediscrepancy.r $filename $tag
