#!/usr/bin/env python2

import numpy as np
from sys import argv
from periodic import element

coords = open(argv[1],"r")
coords_lines = coords.read().splitlines()
natoms = int(coords_lines[0])

symbols=[]
nelectrons=0
for line in coords_lines[2:]:
    nelectrons+= element(line.split()[0]).atomic
    symbols.append(line.split()[0])


print "There are {} atoms and {} electrons".format(natoms,nelectrons)

for atom in np.unique(symbols):
    print "{}: {} atoms x {} e".format(atom,symbols.count(atom),element(atom).atomic)
