#!/bin/bash

filename=$1
title=$2

vcf2tsv \
    | tsvsplit \
        QUAL \
        AC \
        AF \
        TS \
    | tf2binary \
    | vcfplottstv.r $filename $title
