% -*- coding: utf-8 -*-

# VCFLib API

This document describes the VCFLIB API as it is used by the vcflib modules and the [python ffi](./pyvcflib.md).

VCFLIB contains a lot of functionality, but the basis of going through a VCF file and fetching record (by record) information is visible in all modules. A recent example can be found in [vcfwave](../../src/vcfwave.cpp).

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
  cout << var.name << " " << var.pos << endl;
}
```

## Other fields

## Output a VCF record
