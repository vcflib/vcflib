% VCFFIXUP(1) vcffixup (vcflib) | vcffixup (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcffixup**

# SYNOPSIS

**vcffixup** <vcf file>

# DESCRIPTION

Generates a VCF stream where AC and NS have been generated for each record using sample genotypes



# OPTIONS

```




Count the allele frequencies across alleles present in each record in the VCF file. (Similar to vcftools --freq.)

Uses genotypes from the VCF file to correct AC (alternate allele count), AF
(alternate allele frequency), NS (number of called), in the VCF records.  For
example:

    % vcfkeepsamples file.vcf NA12878 | **vcffixup** - | vcffilter -f "AC > 0"

Would downsample file.vcf to only NA12878, removing sites for which the sample
was not called as polymorphic.

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

[vcffixup.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcffixup.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
