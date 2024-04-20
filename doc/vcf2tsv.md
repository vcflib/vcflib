% VCF2TSV(1) vcf2tsv (vcflib) | vcf2tsv (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcf2tsv**

# SYNOPSIS

**vcf2tsv** [-n null_string] [-g] [vcf file]

# DESCRIPTION

Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table. Specifying -g will output one line per sample with genotype information. When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index



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

[vcf2tsv.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2tsv.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
