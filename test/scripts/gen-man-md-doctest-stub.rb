#!/usr/bin/env ruby
#
# Internal helper script to create a markdown file that can be used
# for testing and generating man pages using and erb template
#
# by Pjotr Prins (C) 2020

require 'erb'
require 'date'

VERSION=`cat ../VERSION`.strip

cmd = ARGV.shift
$stderr.print("--- Parsing the output of #{cmd} #{VERSION}\n")

out = `../build/#{cmd} -h 2>&1 `
$stderr.print out

$stderr.print("\n=========================================\n")

lines = out.split(/\n/)
l1 = lines.shift
descr = l1.gsub(/#{VERSION}\s+/,"")
descr = descr.sub(/vcflib/,"VCF")
l2 = lines.shift
l2.length ==0 or raise "Description <#{l2}> should be one line"

l3 = lines.shift
l3 = l3.sub(/..\/build\//,"")
usage_full = l3
usage = l3.sub(/usage:\s+/,"")

lines = lines.map { |x| x == "options:" ? "" : x }

options  = lines.join("\n")

d = DateTime.now

year = d.year

renderer = ERB.new(File.read('scripts/template.erb'))
print output = renderer.result()
