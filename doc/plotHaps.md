% PLOTHAPS(1) plotHaps (vcflib) | plotHaps (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**plotHaps**

# SYNOPSIS



# DESCRIPTION

**plotHaps** provides the formatted output that can be used with 'bin/plotHaplotypes.R'.



# OPTIONS

```


Output : haplotype matrix and positions

**plotHaps**  --target 0,1,2,3,4,5,6,7  --file my.phased.vcf.gz                                                           

required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF column s        
required: r,region     -- argument: a tabix compliant genomic range : "seqid:start-end" or "seqid"                          
required: f,file       -- argument: proper formatted phased VCF file                                                            
required: y,type       -- argument: genotype likelihood format: PL,GP,GP                                                        

Type: statistics

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

[plotHaps.cpp](https://github.com/vcflib/vcflib/blob/master/src/plotHaps.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
