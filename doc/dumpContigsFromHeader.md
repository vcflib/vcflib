% DUMPCONTIGSFROMHEADER(1) dumpContigsFromHeader (vcflib) | dumpContigsFromHeader (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**dumpContigsFromHeader**

# SYNOPSIS

**dumpContigsFromHeader** file

# DESCRIPTION

Dump contigs from header





# EXAMPLES

```

Example:

    **dumpContigsFromHeader** samples/scaffold612.vcf

    ##contig=<ID=scaffold4,length=1524>
    ##contig=<ID=scaffold12,length=56895>
    (...)

    output

    scaffold4       1524
    scaffold12      56895
    (...)

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

[dumpContigsFromHeader.cpp](https://github.com/vcflib/vcflib/blob/master/src/dumpContigsFromHeader.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
