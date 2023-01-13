% VCFAFPATH(1) vcfafpath (vcflib) | vcfafpath (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfafpath**

# SYNOPSIS

**vcfafpath** <vcf file>

# DESCRIPTION

Display genotype paths





# EXAMPLES

```

Example:

    **vcfafpath** samples/scaffold612.vcf

```

T -> A
A -> G
T -> C
C -> A
C -> T
A -> G
T -> C
G -> C
C -> CAGA
A -> G
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

[vcfafpath.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfafpath.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
