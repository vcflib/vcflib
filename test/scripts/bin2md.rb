#!/usr/bin/env ruby
#
# Internal helper script to create markdown docs and man files from
# the --help command of those binaries
#
# The --index option creates the index page
#
#    bin2md [--index] erbtemplate [binary]
#
# by Pjotr Prins (C) 2020-2023
#
# The rules are:
#
#    usage: can be multiline block anywhere in the output
#
# After removing that the next block is the DESCRIPTION.
# Onle line create a TYPE (see TYPES).
# The rest are the OPTIONS.

require 'erb'
require 'date'
require 'open3'

# Section headers for the main index page
TYPES = ["filter","metrics","phenotype","genotype","transformation","statistics"]

=begin
Example of such a section in C++:

if (argc == 2) {
  string h_flag = argv[1];
  if (h_flag == "-h" || h_flag == "--help") {
      cerr << R"(
Generate a random VCF file

Usage: vcfrandom

Example:

    vcfrandom


Type: statistics

      )";
      exit(1);
    }
  }

=end


VERSION=`cat ./VERSION`.strip
template = ARGV.shift
is_man = false # creating man pages?
if template == "--man"
  is_man = true
  template = ARGV.shift
end
create_index = false
if template == "--index"
  create_index = true
  template = ARGV.shift
  index = []
end
search = ARGV.shift

bindir = './build'
$stderr.print("--- Parsing the bin files in #{bindir} for #{VERSION}\n")

d = DateTime.now
year = d.year

# This code walks every binary that was generated, runs it and parses the
# output to create a markdown file.
Dir.glob(bindir+'/*').sort.each do |bin|
  if !File.directory?(bin) and File.executable?(bin)
    if search and bin !~ /#{search}/
      next
    end
    next if not File.executable?(bin) or bin =~ /\.(so|a)$/
    cmd = File.basename(bin)
    help_cmd = cmd + " -h"
    $stderr.print("    "+bin+"\n")
    stdout, stderr, status = Open3.capture3(bin+" -h")
    out = stderr + stdout
    if out == ""
      help_cmd = cmd
      stdout, stderr, status = Open3.capture3(bin)
      out = stderr + stdout
    end
    # out = ,:encoding => 'UTF-8'
    out = out.encode('UTF-8')
    # $stderr.print(out)
    lines = (out).split("\n")
    lines = lines.map{|l| l.gsub(/#{Regexp.escape(cmd)}/,"**#{cmd}**")}
    lines = lines.map{|l| l.gsub(/\.+\/build\//,"")}
    lines = lines.map{|l| l.gsub(/INFO: help:?/,"")}
    lines = lines.map{|l| l.gsub(/INFO: description:/,"")}
    lines = lines.map{|l| l.gsub(/INFO:\s+/,"")}
    pydoc_full = lines.map{|l| l=="" ? '>' : l }.join("\n")
    in_usage = false
    has_options = nil
    has_example = false
    usage = []
    other = []
    example = []
    lines.shift while lines[0] == ""
    lines.each do | l |
      break if l == "------------------------------------------------------"
      if l =~ /usage/i
        in_usage = true
      end
      if in_usage
        if l == "" or l =~ /^Output/i
          in_usage = false
          next
        end
        usage << l
      else
        if l =~ /^Example:/
          has_example = true
        end
        if has_example
          example << l
        else
          other << l
        end
      end
    end
    descr = []
    rest = other
    type = "unknown"
    (other+example).each do | l |
      if l =~ /type: (\S+)/i
        type = $1
        raise "Unknown type #{type} for #{cmd}" if !TYPES.include?(type)
        break
      end
    end

    other.each do | l |
      break if l == "" or l =~ /^Output/i or l =~ /Options/i
      descr << l
    end

    if descr == []
      lineno = 0
      rr = rest.reverse
      rr.each_with_index do | l,i |
        if l != "" and l !~ /^Type/i
          lineno = i
          break
        end
      end
      rr = rr[lineno..-1]
      rr.each do | l |
        if descr.length and l == "" or l =~ /^\s/
          descr = descr.reverse
          break
        end
        descr << l
      end
      rest = rr.drop(descr.size).reverse
    else
      rest = rest.drop(descr.size)
    end

    body = rest.join("\n")
    has_options = true if body != ""
    usage = usage.join(" ").gsub(/#{VERSION}\s+/,"")
    usage = usage.gsub(/usage:\s+/i,"")
    usage = usage.gsub(/\s+/," ").strip
    pydoc_full = pydoc_full.gsub(/#{VERSION}\s+/,"")
    pydoc_full = pydoc_full.gsub(/\.\.\/build\//,"")
    # pydoc_full = pydoc_full.gsub(/vcflib/,"VCF")
    descr = descr.join(" ").gsub(/#{VERSION}\s+/,"")
    descr = descr.sub(/vcflib/,"VCF")
    descr = descr.gsub(/\s+/," ").strip
    example = example.join("\n")
    # print("HELP:",help_cmd,"\n")
    # print("DESCRIPTION:",descr,"\n")
    # print("USAGE:",usage,"\n")
    # print("TYPE:",type,"\n")
    # print("BODY:",body,"\n")

    if create_index
      rec = {
        cmd: cmd,
        type: type,
        descr: descr
      }
      index << rec
    else
      b = binding
      renderer = ERB.new(File.read(template))

      File.open("./doc/#{cmd}.md","w") { |f|
        f.print renderer.result(b)
      }
    end
  end
end

if create_index
  require 'ostruct'

  renderer = ERB.new(File.read("./test/scripts/index-item.erb"))
  File.open("./doc/vcflib.md","w") { |f|
    f.print <<HEADER
% vcflib(1) vcflib | vcflib (index)
% Erik Garrison and vcflib contributors

# NAME

**vcflib** index

# DESCRIPTION

vcflib contains tools and libraries for dealing with the Variant Call
Format (VCF) which is a flat-file, tab-delimited textual format
intended to describe reference-indexed variations between
individuals.

VCF provides a common interchange format for the description of
variation in individuals and populations of samples, and has become
the defacto standard reporting format for a wide array of genomic
variant detectors.

vcflib provides methods to manipulate and interpret sequence variation
as it can be described by VCF. It is both:

* an API for parsing and operating on records of genomic variation as it can be described by the VCF format,
* and a collection of command-line utilities for executing complex manipulations on VCF files.

The API itself provides a quick and extremely permissive method to
read and write VCF files. Extensions and applications of the library
provided in the included utilities (*.cpp) comprise the vast bulk of
the library's utility for most users.

<!--
  Created with ./scripts/bin2md.rb --index
-->

HEADER
    TYPES.each do | type |
      f.print   %{
## #{type}

| #{type} command | description |
| :-------------- | :---------- |
}

      index.each do | rec |
        rec = OpenStruct.new(rec)
        if rec.type == type
          b = binding
          f.print renderer.result(b)
        end
      end
    end
    github = "https://github.com/vcflib/vcflib"
    f.print <<FOOTER

# SOURCE CODE

See the source code repository at #{github}

# CREDIT

Citations are the bread and butter of Science.  If you are using this
software in your research and want to support our future work, please
cite the following publication:

Please cite:

[A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009123).
Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P (2022), PLoS Comput Biol 18(5): e1009123. https://doi.org/10.1371/journal.pcbi.1009123


# LICENSE

Copyright 2011-#{year} (C) Erik Garrison and vcflib contributors. MIT licensed.

FOOTER
  }
end
