#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys

coords = open(argv[1],"r")
outfile = open("stacked_" + argv[1],"w")
coords_lines = coords.read().splitlines()
natoms = int(coords_lines[0])
symbol,x,y,z = [],[],[],[]
print "Define the placement of the second monomer with respect to the first, in Angstroms\n"
xshift = float(raw_input("Shift on X-axis: "))
yshift = float(raw_input("Shift on Y-axis: "))
zshift = float(raw_input("Shift on Z-axis: "))

for line in coords_lines[2:natoms+2]:
    coord = line.split()
    symbol.append(coord[0])
    x.append(float(coord[1]))
    y.append(float(coord[2]))
    z.append(float(coord[3]))

outfile.write("{0}\n\n".format(natoms*2))
for i in range(0,natoms):
    outfile.write("{0} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(symbol[i],x[i],y[i],z[i]))
for i in range(0,natoms):
    outfile.write("{0} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(symbol[i],x[i]+xshift,y[i]+yshift,z[i]+zshift))    
        