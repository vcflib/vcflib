#!/usr/bin/env ruby
#
# Internal helper script to create a markdown file from scripts
#
# by Pjotr Prins (C) 2020

require 'erb'
require 'date'
require 'open3'

# cerr << endl << "Type: statistics" << endl << endl;
# cerr << endl << "Type: transformation" << endl << endl;
TYPES = ["filter","metrics","phenotype","genotype","transformation","statistics"]

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

bindir = './scripts'
$stderr.print("--- Parsing the script files in #{bindir} for #{VERSION}\n")

d = DateTime.now
year = d.year

print   %{
| script | description |
| :-------------- | :---------- |}

Dir.glob(bindir+'/*').each do|bin|
  if !File.directory?(bin) and File.executable?(bin)
    if search and bin !~ /#{search}/
      next
    end
    cmd = File.basename(bin)
    # $stderr.print(cmd,"\n")
    # Get header
    descr = []
    File.open(bin) do | f |
      f.each_line { | line |
        next if line =~ /^\#!/
        break if line !~ /^\#/ or line.strip == "#"
        line = line[0..2].downcase + line[3..-1]
        descr << line
      }
    end
    descr = descr.map{|l| l.sub(/^#\s+/,"").strip }
    print   %{
| [#{cmd}](./scripts/#{cmd}) | #{descr.join(" \\\n").strip} |}
  end
end
print("\n")
