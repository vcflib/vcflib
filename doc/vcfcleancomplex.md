% VCFCLEANCOMPLEX(1) vcfcleancomplex (vcflib) | vcfcleancomplex (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcleancomplex**

# SYNOPSIS

**vcfcleancomplex** <vcf file>

# DESCRIPTION

Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change.



# OPTIONS

```


Generate a VCF stream in which 'long' non-complexalleles have their position corrected.
assumes that VCF records can't overlap 5'->3'

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

[vcfcleancomplex.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcleancomplex.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
