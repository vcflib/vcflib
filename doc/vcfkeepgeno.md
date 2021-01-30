% VCFKEEPGENO(1) vcfkeepgeno (vcflib) | vcfkeepgeno (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfkeepgeno**

# SYNOPSIS

**vcfkeepgeno** <vcf file> [FIELD1] [FIELD2] ...

# DESCRIPTION

Reduce file size by removing FORMAT fields not listed on the command line from sample specifications in the output



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

[vcfkeepgeno.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfkeepgeno.cpp)

# LICENSE

Copyright 2011-2021 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
