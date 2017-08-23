#!/usr/bin/env python

from sys import argv

# This code takes an xyz file containing multiple geometries and splits them
# into individual xyz files. Each geometry must have the samme number of atoms
# as the initial structure

geomfile=open(argv[1],"r").read().splitlines()
natoms=int(geomfile[0])
nlines=natoms+2
count=1
for geom in range(len(geomfile)/(nlines)):
    outfile=open(argv[1]+"-"+str(count)+".xyz","w")
    for line in geomfile[0+(geom*nlines):nlines+(geom*nlines)]:
        outfile.write("{}\n".format(line))
    count+=1
