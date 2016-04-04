#!/bin/bash

tag=$1

vcf2tsv \
    | tsvsplit \
        QUAL \
        AC \
        $tag.has_variant \
        $tag.site.alternate_negative_discrepancy \
        $tag.site.alternate_positive_discrepancy \
        $tag.genotypes.alternate_count \
        $tag.site.non_reference_sensitivity.count \
        $tag.site.non_reference_sensitivity.normalizer \
        $tag.site.non_reference_discrepancy.count \
        $tag.site.non_reference_discrepancy.normalizer \
    | tf2binary \
    | vcfprintaltdiscrepancy.r $tag
