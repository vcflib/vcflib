% VCFANNOTATEGENOTYPES(1) vcfannotategenotypes (vcflib) | vcfannotategenotypes (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfannotategenotypes**

# SYNOPSIS

**vcfannotategenotypes** <annotation-tag> <vcf file> <vcf file>

# DESCRIPTION

Examine genotype correspondence. Annotate genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant.



# OPTIONS

```


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

[vcfannotategenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotategenotypes.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
