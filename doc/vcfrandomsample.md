% VCFRANDOMSAMPLE(1) vcfrandomsample (vcflib) | vcfrandomsample (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfrandomsample**

# SYNOPSIS

**vcfrandomsample** [options] [<vcf file>]

# DESCRIPTION

Randomly sample sites from an input VCF file, which may be provided as stdin. Scale the sampling probability by the field specified in KEY. This may be used to provide uniform sampling across allele frequencies, for instance.



# OPTIONS

```

options:
    -r, --rate RATE          base sampling probability per locus
    -s, --scale-by KEY       scale sampling likelihood by this Float info field
    -p, --random-seed N      use this random seed (by default read from /dev/random)
    -q, --pseudorandom-seed  use a pseudorandom seed (by default read from /dev/random)


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

[vcfrandomsample.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfrandomsample.cpp)

# LICENSE

Copyright 2011-2021 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
