# vcflib performance

## vcfwave

The early edition of WFA2-lib did not automatically compile for parallel. In the 'build' dir:

```
/usr/bin/time -v ./vcfwave ../samples/grch38#chr8_36353854-36453166-bcftools-normalised.vcf > /dev/null
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:09.18
```

```
make BUILD_WFA_PARALLEL=1 BUILD_TOOLS=0  BUILD_TOOLS=0 BUILD_EXAMPLES=0 CC=gcc CC_FLAGS=-fPIC setup asan lib_wfa
```

I made sure to build vcflib with OPENMP=ON and libgomp.so.1 is linked. Also WFA builds with openmp:

```
CMakeFiles/wfa-EXT.dir/build.make:      cd /export/local/home/wrk/iwrk/opensource/code/pangenome/vcflib/contrib/WFA2-lib && $(MAKE) clean BUILD_WFA_PARALLEL=1 BUILD_TOOLS=0 BUILD_EXAMPLES=0 CC=gcc CC_FLAGS=-fPIC setup lib_wfa
```

This happens when running `cmake -DOPENMP=ON ..`, but, still, no multi-core `vcfwave`!

When setting

```C
wavefront_aligner_set_max_num_threads(wf_aligner,4);
```

We got a speedup:

        Percent of CPU this job got: 242%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.24

### Big files

This is a 30Mb file after compression:

```
tux02:~/tmp/vcflib/build$ /usr/bin/time -v ./vcfwave ../samples/chr18.grch38.vcfbub.a100k.vcf.gz > test.vcf
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:17:16
```

Trying 16 threads

```
tux02:~/tmp/vcflib/build$ /usr/bin/time -v ./vcfwave ../samples/chr18.grch38.vcfbub.a100k.vcf.gz -t 16 > test16.vcf
        Percent of CPU this job got: 464%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:16:30
```

hmm. Not so exciting. Looks like we need to parallelize at a different level. Memory used was about 6Gb, and that is not bad.

The good news is that the output is identical.

## vcfcreatemulti

This program is actually slower than vcfwave.
