% VCFEVENREGIONS(1) vcfevenregions (vcflib) | vcfevenregions (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfevenregions

# SYNOPSIS

usage: ./build/vcfevenregions [options] <vcf file>

# DESCRIPTION

options: -f, --fasta-reference REF FASTA reference file to use to obtain primer sequences. -n, --number-of-regions N The number of desired regions. -p, --number-of-positions N The number of positions per region. -o, --offset N Add an offset to region positioning, to avoid boundary related artifacts in downstream processing. -l, --overlap N The number of sites to overlap between regions. Default 0. -s, --separator SEQ Specify string to use to separate region output. Default '-'

# OPTIONS

```


Generates a list of regions, e.g. chr20:10..30 using the variant
density information provided in the VCF file to ensure that the regions have
even numbers of variants.  This can be use to reduce the variance in runtime
when dividing variant detection or genotyping by genomic coordinates.

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfevenregions.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfevenregions.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
