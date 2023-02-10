// Python ffi calls C++ functions
//
// Copyright © 2022-2023 Pjotr Prins

#include "Variant.h"
#include "vcf-wfa.h"

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace vcflib;

PYBIND11_MODULE(pyvcflib, m)
{
  m.doc() = "This is a Python binding of C++ vcflib Library with WFA";

  // Wavefront

  py::class_<affine2p_penalties_t>(m, "affine2p_penalties_t", "WFA affine2 settings")
      .def_readwrite("match", &affine2p_penalties_t::match)
      .def_readwrite("mismatch", &affine2p_penalties_t::mismatch)
      .def_readwrite("gap_opening1", &affine2p_penalties_t::gap_opening1)
      .def_readwrite("gap_extension1", &affine2p_penalties_t::gap_extension1)
      .def_readwrite("gap_opening2", &affine2p_penalties_t::gap_opening2)
      .def_readwrite("gap_extension2", &affine2p_penalties_t::gap_extension2)
      ;
  py::class_<wavefront_aligner_attr_t>(m, "wavefront_aligner_attr_t", "WFA settings")
      .def_readwrite("distance_metric", &wavefront_aligner_attr_t::distance_metric)
      .def_readwrite("affine2p_penalties", &wavefront_aligner_attr_t::affine2p_penalties)
      .def_readwrite("alignment_scope", &wavefront_aligner_attr_t::alignment_scope)
      ;
  m.attr("wavefront_aligner_attr_default") = wavefront_aligner_attr_default;
  py::enum_<distance_metric_t>(m, "distance_meric_t")
      .value("gap_affine_2p", gap_affine_2p)
      ;
  py::enum_<alignment_scope_t>(m, "alignment_scope_t")
      .value("compute_alignment", compute_alignment)
      ;


  // Main VCFlib
  py::class_<VariantAllele>(m, "VariantAllele", "VCF alleles")
      .def_readonly("position", &VariantAllele::position)
      .def_readonly("ref", &VariantAllele::ref)
      .def_readonly("alt", &VariantAllele::alt)
      ;
  py::class_<Variant>(m, "Variant", "VCF record")
      .def(py::init<VariantCallFile &>() )
      .def_readwrite("name", &Variant::sequenceName)
      .def_readwrite("pos", &Variant::position)
      .def_readwrite("ref", &Variant::ref)
      .def_readwrite("alt", &Variant::alt)
      .def_readwrite("alleles", &Variant::alleles)
      .def_readonly("info", &Variant::info)
      .def_readonly("sampleNames", &Variant::sampleNames)
      .def_readonly("samples", &Variant::samples)
      .def("legacy_parsedAlternates", &Variant::legacy_parsedAlternates)
      ;

  py::class_<WfaVariant, Variant>(m, "WfaVariant", "WFA VCF record")
      .def(py::init<VariantCallFile &>() )
      .def("parsedAlternates", &WfaVariant::parsedAlternates);

  py::class_<VariantCallFile>(m, "VariantCallFile", "VCF file")
      .def(py::init())
      .def("openFile",&VariantCallFile::openFile,"Open the VCF")
      .def("getNextVariant",&VariantCallFile::getNextVariant,"Iterate VCF records")
      ;
}
