#!/usr/bin/env python

import sys
import re
import sqlite3

if len(sys.argv) < 2:
    print "usage", sys.argv[0], " [dbname]"
    print "reads VCF on stdin, and writes output to a sqlite3 db [dbname]"
    exit(1)

dbname = sys.argv[1]

# parse the header
# into a mapping from tag -> type

infotypes = {}
infonumbers = {}

for line in sys.stdin:
    if line.startswith('##INFO'):
        #<ID=XRS,Number=1,Type=Float,
        i = re.search("ID=(.*?),", line)
        n = re.search("Number=(.*?),", line)
        t = re.search("Type=(.*?),", line)
        if i and n and t:
            id = i.groups()[0]
            number = n.groups()[0]
            if number == "A":
                number = -1
            elif number == "G" or int(number) > 1:
                # unclear how to deal with these
                continue
            else:
                number = int(number)
            typestr = t.groups()[0]
            infotypes[id] = typestr
            infonumbers[id] = number
        else:
            continue
    elif line.startswith('##'):
        continue
    else:
        break # header line, sample names etc.

# write the table schema

infotype_to_sqltype = {}
infotype_to_sqltype["Flag"] = "boolean"
infotype_to_sqltype["Integer"] = "integer"
infotype_to_sqltype["Float"] = "real"
infotype_to_sqltype["String"] = "text"

tablecmd = """create table alleles"""
specs = ["CHROM text",
        "POS integer",
        "ID text",
        "REF text",
        "ALT text",
        "QUAL real",
        "FILTER text"]

sorted_fields = sorted(infotypes.keys())
for field in sorted_fields:
    infotype = infotypes[field]
    sqltype = infotype_to_sqltype[infotype]
    field = field.replace(".", "_") # escape periods, which are not allowed
    specs.append(field + " " + sqltype)

tablecmd += " (" + ", ".join(specs) + ")"

conn = sqlite3.connect(dbname)
conn.execute(tablecmd)

# for each record
# parse the record
# for each allele

for line in sys.stdin:
    fields = line.split('\t')
    chrom, pos, id, ref, alts, qual, filter, info = fields[:8]
    alts = alts.split(",")
    altindex = 0
    chrom = "\'" + chrom + "\'"
    id = "\'" + id + "\'"
    ref = "\'" + ref + "\'"
    filter = "\'" + filter + "\'"
    for alt in alts:
        alt = "\'" + alt + "\'"
        info_values = {}
        for pair in info.split(";"):
            if pair.find("=") is not -1:
                pair = pair.split("=")
                key = pair[0]
                value = pair[1]
                if not infonumbers.has_key(key):
                    continue
                if infonumbers[key] == -1:
                    values = value.split(",")
                    value = values[altindex]
                info_values[key] = value
            else:
               # boolean flag
                info_values[pair] = "1"
        ordered_insertion = []
        for field in sorted_fields:
            value = "null"
            if info_values.has_key(field):
                value = info_values[field]
                if infotypes[field] == "String":
                    value = "\'" + value + "\'"
            else:
                # missing flag means "false" for that flag
                if infotypes[field] == "Flag":
                    value = "0"
            ordered_insertion.append(value)
        cmd = "insert into alleles values (" \
            + ", ".join([chrom, pos, id, ref, alt, qual, filter]) \
            + ", " \
            + ", ".join(ordered_insertion) + ")"
        conn.execute(cmd)
        altindex += 1

conn.commit()

# TODO ignoring samples (for now)

# add indexes everywhere?

conn.close()
