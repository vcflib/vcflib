% VCFUNIQ(1) vcfuniq (vcflib) | vcfuniq (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcfuniq**

# SYNOPSIS

**vcfuniq** <vcf file>

# DESCRIPTION

List unique genotypes. Similar to GNU uniq, but aimed at VCF records. **vcfuniq** removes records which have the same position, ref, and alt as the previous record on a sorted VCF file. Note that it does not adjust/combine genotypes in the output, but simply takes the first record. See also vcfcreatemulti for combining records.



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

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
