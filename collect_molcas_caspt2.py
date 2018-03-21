#!/usr/bin/env python

import sys

with open(sys.argv[1],"r") as infile:
    energies=[]
    for line in infile:
        if "MS-CASPT2 Root" in line:
            energies.append(line.split()[5])
if not energies:
    sys.exit("Couldn't find any MS-CASPT2 roots!")
else:
    energies_ev=[(i-energies[0])*27.2114 for i in energies]

print("Energies of MS-CASPT2 states:\n")
for i,energy in enumerate(energies_ev):
    print("{0}: {1:.2f} eV".format(i,energy))
