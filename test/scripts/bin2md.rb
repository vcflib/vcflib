#!/usr/bin/env ruby
#
# Internal helper script to create a markdown file from binaries
#
# The --index option creates the index page
#
#    bin2md [--index] erbtemplate [binary]
#
# by Pjotr Prins (C) 2020
#
# The rules are simple USAGE:
#
#    usage: can be multiline block anywhere in the output
#
# After removing that the next block is the DESCRIPTION.
# Onle line create a TYPE (see TYPES).
# The rest are the OPTIONS.

require 'erb'
require 'date'
require 'open3'

# Add with cerr << endl << "Type: statistics" << endl << endl;
# Add with cerr << endl << "Type: transformation" << endl << endl;
TYPES = ["filter","transformation","statistics","metrics","phenotype","genotype"]

=begin
  if (argc == 2) {
    string h_flag = argv[1];

    if (argc == 2 && (h_flag == "-h" || h_flag == "--help")) {
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
create_index = false
if template == "--index"
  create_index = true
  template = ARGV.shift
end
search = ARGV.shift

bindir = './build'
$stderr.print("--- Parsing the bin files in #{bindir} for #{VERSION}\n")

Dir.glob(bindir+'/*').each do|bin|
  if !File.directory?(bin) and File.executable?(bin)
    if search and bin !~ /#{search}/
      next
    end
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

    body = rest.drop(descr.size).join("\n")
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

    d = DateTime.now

    year = d.year

    b = binding
    renderer = ERB.new(File.read(template))

    File.open("./doc/#{cmd}.md","w") { |f|
      f.print renderer.result(b)
    }

  end
end
