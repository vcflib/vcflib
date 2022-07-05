% -*- coding: utf-8 -*-

# VCFLib API

This document describes the VCFLIB API as it is used by the vcflib modules and the [python ffi](./pyvcflib.md).

VCFLIB contains a lot of functionality, but the basis of going through a VCF file and fetching record (by record) information is straightforward and visible in all modules. A recent example can be found in [vcfwave](../../src/vcfwave.cpp).

## Open the VCF file

```C++
VariantCallFile variantFile;

if (optind < argc) {
    string filename = argv[optind];
    variantFile.open(filename);
} else {
    variantFile.open(std::cin);
}

if (!variantFile.is_open()) {
    return 1;
}
```

## Read records

The following will parse the records and you can print out the first two fields with

```C++
Variant var(variantFile);
while (variantFile.getNextVariant(var)) {
  cout << var.sequenceName << " " << var.position << endl;
}
```

## Other fields

In the file [Variant.h](https://github.com/vcflib/vcflib/blob/master/src/Variant.h) the Variant class is defined with fields/accessors, such as

```C++
    string sequenceName;
    long position;
    long zeroBasedPosition(void) const;
    string id;
    string ref;
    vector<string> alt;      // a list of all the alternate alleles present at this locus
    vector<string> alleles;  // a list all alleles (ref + alt) at this locus
```

See above read records example.

## Output a VCF record

The default string outputter of the Variant class outputs a VCF record using the field that are defined:

```C++
Variant var(variantFile);
while (variantFile.getNextVariant(var)) {
  cout << var << endl;
}
```

will output VCF.

## Mirroring in the Python FFI

The Python FFI follows this API though some accessors may be renamed. See
[pythonffi.cpp](pythonffi.cpp) and [pyvcflib.md](pyvcflib.md).
