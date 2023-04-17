% VCFREMAP(1) vcfremap (vcflib) | vcfremap (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfremap**

# SYNOPSIS

**vcfremap** [options] [<vcf file>]

# DESCRIPTION

For each alternate allele, attempt to realign against the reference with lowered gap open penalty. If realignment is possible, adjust the cigar and reference/alternate alleles. Observe how different alignment parameters, including context and entropy-dependent ones, influence variant classification and interpretation.



# OPTIONS

```

options:
    -w, --ref-window-size N      align using this many bases flanking each side of the reference allele
    -s, --alt-window-size N      align using this many flanking bases from the reference around each alternate allele
    -r, --reference FILE         FASTA reference file, required with -i and -u
    -m, --match-score N          match score for SW algorithm
    -x, --mismatch-score N       mismatch score for SW algorithm
    -o, --gap-open-penalty N     gap open penalty for SW algorithm
    -e, --gap-extend-penalty N   gap extension penalty for SW algorithm
    -z, --entropy-gap-open       use entropy scaling for the gap open penalty
    -R, --repeat-gap-extend N    penalize non-repeat-unit gaps in repeat sequence
    -a, --adjust-vcf TAG         supply a new cigar as TAG in the output VCF


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfremap.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfremap.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
