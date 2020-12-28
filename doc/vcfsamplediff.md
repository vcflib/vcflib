% VCFSAMPLEDIFF(1) vcfsamplediff (vcflib) | vcfsamplediff (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfsamplediff**

# SYNOPSIS

**vcfsamplediff** [options] <tag> <sample> <sample> [ <sample> ... ] <vcf file>

# DESCRIPTION

Establish putative somatic variants using reported differences between germline and somatic samples. Tags each record where the listed sample genotypes differ with <tag>. The first sample is assumed to be germline, the second somatic. Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples.



# OPTIONS

```


options:
    -s --strict     Require that no observations in the germline support the somatic alternate.


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

[vcfsamplediff.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsamplediff.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
