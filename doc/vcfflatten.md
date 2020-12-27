% VCFFLATTEN(1) vcfflatten (vcflib) | vcfflatten (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfflatten**

# SYNOPSIS

**vcfflatten** [file]

# DESCRIPTION

Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin.





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfflatten.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfflatten.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
