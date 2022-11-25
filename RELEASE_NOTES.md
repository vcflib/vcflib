For contributions see
[contributors](https://github.com/vcflib/vcflib/graphs/contributors)
and
[commits](https://github.com/vcflib/vcflib/commits/master).

## TODO

- [ ] vcfcreatemulti merge multiple rows
      + [X] rewrite vcfcreatemulti using zig
      + [.] merge genotypes correctly, with tests
      + [ ] adjust info and genotypes for variants that have multiple alts already (now errors)
      + [X] handle phase
      + [ ] check for indels which are really the same
      + [ ] combine vcfwave duplicated functionality
- [ ] RELEASE 1.0.5
- [ ] vcfuniq combine genotypes from non-unique records
- [ ] vcfwave check polyploids?
- [ ] RELEASE 1.0.6

## ChangeLog v1.0.4 (20220829)

This is a release with some major changes.

Important changes:

- vcfwave is introduced and vcfallelicprimitimes is now considered obsolete
- INFO fields output order is now the same as on input parsing
- vcfwave check merging of genotypes - write tests
- vcfwave recompute AC, AFs from merged record
- introduced the zig compiler with vcfcreatemulti.cpp as a first target (use cmake ZIG=OFF to disable). At this point the zig version (-n switch) gives identical results.

Introduction of O(n) wavefront algorithm WF to replace O(n^2) Smith-Waterman SW. Note that the output is different from the original SW implementation. SW is still optionally available but considered obsolete. Use the bi-directional vcfwave instead of vcfallelicprimitives.

+ Added realignment using the wavefront algorithm (now the default). See [vcfwave](./doc/vcfwave.md) (thank you Erik Garrison https://github.com/ekg and Santiago Marco-Sola  https://github.com/smarco).
+ Support longer read inversions in vcfwave!
+ vcfallelicprimitives now considered legacy/obsolete
+ Improved CMake configuration
+ vcflib compiles with both gcc and clang++ and tests pass, see [guix-clang.scm](./guix-clang.scm) - mind that git submodules such as WFA2-lib still override to gcc
+ Fixed local build for tabixpp+htslib - note that htslib should be an upstream released version (currently 1.15.1). Unfortunately git submodule does not handle tags.
+ Fix for -L switch for vcfallelicprimitives
+ Added libasan and lto support
+ Removed useless googletest submodule
+ Added python bindings with pybind11
+ Added python testing framework
+ Added tabixpp back in as a submodule, fixes https://github.com/vcflib/vcflib/issues/305
+ Optimizations and bug fixes. (thanks @mphschmitt)

## ChangeLog v1.0.3 (20220122)

This is a maintenance release of vcflib.

+ Merge intervaltree changes (thanks @jnunn and @timmassingham)
+ Built with gcc-11
+ Fix issue #251 hapLrt: fix segfault when accessing genotype field. (thanks @mphschmitt)
+ Fix vcfflatten: fix segfault when no 'AF' field is present (#47, thanks @mphschmitt)
+ Fixes on vcfnulldotslashdot #310 (thanks @WinterFor)
+ Fix issue #301: Replace raw pointer usage with std::unique_ptr #306 (thanks @Glebanister)
+ Fix man page installation #321 (thanks @alexreg)
+ Use `guix shell` instead of `guix environment` for development
+ Regenerated online docs
+ README: add matrix badge (removed gitter badge)

## ChangeLog v1.0.2 (20210104)

This is a maintenance release of vcflib, mostly improving the build
system, CI and generating markdown docs as well as man pages.

+ Removed tabixpp and htslib source dependencies, i.e., we are now using
  the distro provided libraries and include files through pkg-config.
  See also the [README](README.md#build-from-source)
+ Removed the tabixpp+htslib git submodules
+ Generalise and document the cmake build system
+ Added tests to the cmake build system and build instructions to README
+ Added support for ARM64 and PowerPC, see #292 (thanks @genisysram and @mr-c)
+ Added github actions for the issue tracker
+ Added githum CI
+ Updated header files in src with copyright/license info, see #16
+ Created markdown [docs](./doc/vcflib.md) and [man pages](./man/) for
  all utilities. Created a script bin2md for markdown generation and
  use pandoc for the man page generation.

## Older changes

For older changes view the git [log](https://github.com/vcflib/vcflib/commits/master).
