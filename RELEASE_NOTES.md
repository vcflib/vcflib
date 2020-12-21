## ChangeLog v1.0.2 (2020????)

This is a maintenance release of Freebayes:

+ Removed tabixpp and htslib source dependencies, i.e., we are now using
  the distro provided libraries and include files through pkg-config.
+ Removed the tabixpp+htslib git submodule
+ Generalise and document the cmake build system
+ Added tests to the cmake build system and build instructions to README
+ Added support for ARM64 and PowerPC, see #292 (thanks @genisysram and @mr-c)
+ Added github actions for the issue tracker
+ Added githum CI
+ Updated header files in src with copyright/license info, see #16

## Older changes

For older changes view the git [log](https://github.com/vcflib/vcflib/commits/master).
