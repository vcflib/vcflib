% VCFROC(1) vcfroc (vcflib) | vcfroc (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfroc

# SYNOPSIS

usage: ./build/vcfroc [options] [<vcf file>]

# DESCRIPTION

options: -t, --truth-vcf FILE use this VCF as ground truth for ROC generation -w, --window-size N compare records up to this many bp away (default 30) -c, --complex directly compare complex alleles, don't parse into primitives -r, --reference FILE FASTA reference file



# OPTIONS

```


Generates a pseudo-ROC curve using sensitivity and specificity estimated against
a putative truth set.  Thresholding is provided by successive QUAL cutoffs.

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfroc.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfroc.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
