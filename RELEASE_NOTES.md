iFor contributions
see
[contributors](https://github.com/vcflib/vcflib/graphs/contributors)
and
[commits](https://github.com/vcflib/vcflib/commits/master).

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
