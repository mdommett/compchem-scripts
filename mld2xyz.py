#!/usr/bin/env python

import sys
infile=open(sys.argv[1],"r").read().splitlines()
outfile=open(sys.argv[2],"w")
for i,line in enumerate(infile):
    if "[GEOMETRIES] (XYZ)" in line:
        start=i+1
    if "[FORCES]" in line:
        stop=i


for line in infile[start:stop]:
    outfile.write("{}\n".format(line))
