% VCFLEFTALIGN(1) vcfleftalign (vcflib) | vcfleftalign (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfleftalign**

# SYNOPSIS

./build/**vcfleftalign** [options] [file]

# DESCRIPTION

options: -r, --reference FILE Use this reference as a basis for realignment. -w, --window N Use a window of this many bp when left aligning (150).



# OPTIONS

```


Left-aligns variants in the specified input file or stdin.  Window size is determined
dynamically according to the entropy of the regions flanking the indel.  These must have
entropy > 1 bit/bp, or be shorter than ~5kb.

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfleftalign.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfleftalign.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
