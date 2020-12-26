% VCFANNOTATEGENOTYPES(1) vcfannotategenotypes (vcflib) | vcfannotategenotypes (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfannotategenotypes

# SYNOPSIS

usage: ./build/vcfannotategenotypes <annotation-tag> <vcf file> <vcf file> annotates genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant.

# DESCRIPTION



# OPTIONS

```



```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfannotategenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotategenotypes.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
