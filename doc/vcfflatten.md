% VCFFLATTEN(1) vcfflatten (vcflib) | vcfflatten (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfflatten**

# SYNOPSIS

**vcfflatten** [options][file] -h --help display this help message and exit. -i --ignore-errors do not flatten locus if 'AF' is not specified.

# DESCRIPTION

Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin.



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

[vcfflatten.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfflatten.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
