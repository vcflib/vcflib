% VCFEVENREGIONS(1) vcfevenregions (vcflib) | vcfevenregions (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfevenregions**

# SYNOPSIS

**vcfevenregions** [options] <vcf file>

# DESCRIPTION

Generates a list of regions, e.g. chr20:10..30 using the variant density information provided in the VCF file to ensure that the regions have even numbers of variants. This can be use to reduce the variance in runtime when dividing variant detection or genotyping by genomic coordinates.



# OPTIONS

```

options:
    -f, --fasta-reference REF    FASTA reference file to use to obtain primer sequences.
    -n, --number-of-regions N    The number of desired regions.
    -p, --number-of-positions N  The number of positions per region.
    -o, --offset N               Add an offset to region positioning, to avoid boundary
                                 related artifacts in downstream processing.
    -l, --overlap N              The number of sites to overlap between regions.  Default 0.
    -s, --separator SEQ          Specify string to use to separate region output.  Default '-'

Type: transformation

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

[vcfevenregions.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfevenregions.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
