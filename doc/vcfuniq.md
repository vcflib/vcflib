% VCFUNIQ(1) vcfuniq (vcflib) | vcfuniq (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcfuniq**

# SYNOPSIS

**vcfuniq** <vcf file>

# DESCRIPTION

List unique genotypes. Like GNU uniq, but for VCF records. Remove records which have the same position, ref, and alt as the previous record.



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

[vcfuniq.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfuniq.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
