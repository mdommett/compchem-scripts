#!/usr/bin/env python

import sys

with open(sys.argv[1],"r") as infile:
    energies=[]
    for line in infile:
        if "Total energy:" in line:
            energies.append(float(line.split()[7]))
if not energies:
    sys.exit("Couldn't find any RASSCF roots!")
else:
    energies_ev=[(i-energies[0])*27.2114 for i in energies]

print("Energies of RASSCF states:")
for i,energy in enumerate(energies_ev):
    print("S{0}: {1:.2f} eV".format(i,energy))
