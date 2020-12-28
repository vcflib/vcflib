% VCFLEFTALIGN(1) vcfleftalign (vcflib) | vcfleftalign (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfleftalign**

# SYNOPSIS

**vcfleftalign** [options] [file]

# DESCRIPTION

Left-align indels and complex variants in the input using a pairwise ref/alt alignment followed by a heuristic, iterative left realignment process that shifts indel representations to their absolute leftmost (5') extent.



# OPTIONS

```


This is the same procedure used in the internal left alignment in
freebayes, and can be used when preparing VCF files for input to
freebayes to decrease positional representation differences between
the input alleles and left-realigned alignments.

options:

        -r, --reference FILE  Use this reference as a basis for realignment.
        -w, --window N        Use a window of this many bp when left aligning (150).

Left-aligns variants in the specified input file or stdin.  Window
size is determined dynamically according to the entropy of the regions
flanking the indel.  These must have entropy > 1 bit/bp, or be shorter
than ~5kb.


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

[vcfleftalign.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfleftalign.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
