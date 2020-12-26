% PLOTHAPS(1) plotHaps (vcflib) | plotHaps (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

plotHaps

# SYNOPSIS



# DESCRIPTION

plotHaps provides the formatted output that can be used with 'bin/plotHaplotypes.R'. Output : haplotype matrix and positions



# OPTIONS

```


plotHaps  --target 0,1,2,3,4,5,6,7  --file my.phased.vcf.gz                                                           

required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF column s        
required: r,region     -- argument: a tabix compliant genomic range : "seqid:start-end" or "seqid"                          
required: f,file       -- argument: proper formatted phased VCF file                                                            
required: y,type       -- argument: genotype likelihood format: PL,GP,GP                                                        

version 1.0.1 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu 

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[plotHaps.cpp](https://github.com/vcflib/vcflib/blob/master/src/plotHaps.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
