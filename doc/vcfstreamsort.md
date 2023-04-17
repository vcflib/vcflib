% VCFSTREAMSORT(1) vcfstreamsort (vcflib) | vcfstreamsort (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfstreamsort**

# SYNOPSIS

**vcfstreamsort** [options] [vcf file]

# DESCRIPTION

Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart.



# OPTIONS

```

options:

    -h, --help             this dialog
    -w, --window N         number of sites to sort (default 10000)
    -a, --all              load all sites and then sort in memory

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

[vcfstreamsort.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfstreamsort.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
