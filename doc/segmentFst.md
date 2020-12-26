% SEGMENTFST(1) segmentFst (vcflib) | segmentFst (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

segmentFst

# SYNOPSIS

usage: segmentFst -s 0.7 -f wcFst.output.txt

# DESCRIPTION

Creates genomic segments (bed file) for regions with high wcFst Output : 8 columns : 1. Seqid 2. Start (zero based) 3. End (zero based) 4. Average Fst 5. Average high Fst (Fst > -s) 6. N Fst values in segment 7. N high fst values in segment 8. Segment length required: -f -- Output from wcFst optional: -s -- High Fst cutoff [0.8]





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[segmentFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/segmentFst.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
