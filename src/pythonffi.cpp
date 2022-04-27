// Python ffi calls C++ functions
//
// Copyright © 2022 Pjotr Prins

#include "Variant.h"

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace vcflib;

PYBIND11_MODULE(pyvcflib, m)
{
  m.doc() = "This is a Python binding of C++ vcflib Library";
  py::class_<Variant>(m, "Variant", "VCF record")
    .def(py::init<VariantCallFile &>() )
    .def_readwrite("name", &Variant::sequenceName)
    .def_readwrite("pos", &Variant::position)
    .def_readwrite("ref", &Variant::ref)
    .def_readwrite("alt", &Variant::alt)
    .def_readwrite("alleles", &Variant::alleles)
    ;
  py::class_<VariantCallFile>(m, "VariantCallFile", "VCF file")
    .def(py::init())
    .def("openFile",&VariantCallFile::openFile,"Open the VCF")
    .def("getNextVariant",&VariantCallFile::getNextVariant,"Iterate VCF records")
    ;
}
