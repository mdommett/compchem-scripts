#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys

molecule = open("molecule","r").read().splitlines()
input_coords = open(argv[1],"r").read().splitlines()
outfile = open(argv[2]+".com","w")


###############
chk = "%chk={0}.chk".format(argv[2])
nproc = "%nproc={0}".format(raw_input("Number of processors: "))
mem = "%mem={0}GB".format(raw_input("Total memory in GB : "))
input_line = "#p ONIOM(wB97XD/6-31G*:amber=(hardfirst))=(EmbedCharge) nosymm opt geom=connectivity"
###############

def group_to_range(group):
    group = ''.join(group.split())
    sign, g = ('-', group[1:]) if group.startswith('-') else ('', group)
    r = g.split('-', 1)
    r[0] = sign + r[0]
    r = sorted(int(__) for __ in r)
    return range(r[0], 1 + r[-1])

def rangeexpand(txt):
    ranges = chain.from_iterable(group_to_range(__) for __ in txt.split(','))
    return sorted(set(ranges))

def read_xyz(txt):
    symbol,x,y,z = [],[],[],[]
    for line in txt[2:]:
        coords = line.split()
	if coords:
		symbol.append(coords[0])
        	x.append(float(coords[1]))
        	y.append(float(coords[2]))
        	z.append(float(coords[3]))
    if len(symbol) == len(x) == len(y) == len(z):
        return symbol,x,y,z
    else:
        print "Error in reading xyz file. Bye..."
        sys.exit

symbol,x,y,z = read_xyz(input_coords)
natoms_inputxyz = len(symbol)
natoms_molecule = int(molecule[0])
Nlevels =  len(molecule)-1

if natoms_inputxyz != natoms_molecule:
    print "There is a problem with the input xyz or the molecule file. The number of atoms don't match."
    sys.exit

if Nlevels == 2:
    H = rangeexpand(molecule[1])
    M = rangeexpand(molecule[2])
    print "{0} atoms in High Level, {1} atoms in Mid level. There {2} atoms in total.".format(len(H),len(M),natoms_molecule)
    if len(H+M) != len(symbol):
        print "There is a problem with the molecule input file. There are {0} in molecule but {1} in the xyz file".format(len(H+M),len(symbol))
        sys.exit

        
if Nlevels == 3:
    H = rangeexpand(molecule[1])
    M = rangeexpand(molecule[2])
    L = rangeexpand(molecule[3])
    print "{0} atoms in High Level, {1} atoms in Mid level {2} atoms in Low level. There {3} atoms in total.".format(len(H),len(M),len(L),natoms_molecule)


print "Writing .com file..."
outfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n".format(chk,nproc,mem,input_line))


for i in H:
    outfile.write("{0:<2}  0 {1:>13.9f} {2:>13.9f} {3:>13.9f} H\n".format(symbol[i-1],x[i-1],y[i-1],z[i-1]))

for i in M:
    outfile.write("{0:<2} 0 {1:>13.9f} {2:>13.9f} {3:>13.9f} M\n".format(symbol[i-1],x[i-1],y[i-1],z[i-1]))

if Nlevels == 3:
    for i in L:
        outfile.write("{0:<2} -1 {1:>13.9f} {2:>13.9f} {3:>13.9f} L\n".format(symbol[i-1],x[i-1],y[i-1],z[i-1]))
    
outfile.write("\n")
