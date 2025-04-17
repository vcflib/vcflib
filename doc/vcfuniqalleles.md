% VCFUNIQALLELES(1) vcfuniqalleles (vcflib) | vcfuniqalleles (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcfuniqalleles**

# SYNOPSIS

**vcfuniqalleles** <vcf file>

# DESCRIPTION

List unique alleles For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files.



# OPTIONS

```


Type: filter

      

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

[vcfuniqalleles.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfuniqalleles.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
