% VCFSTREAMSORT(1) vcfstreamsort (vcflib) | vcfstreamsort (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfstreamsort**

# SYNOPSIS

./build/**vcfstreamsort** [options] [vcf file]

# DESCRIPTION

Sorts the input (either stdin or file) using a streaming sort algorithm. options:



# OPTIONS

```


    -h, --help             this dialog
    -w, --window N         number of sites to sort (default 10000)
    -a, --all              load all sites and then sort in memory

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfstreamsort.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfstreamsort.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
