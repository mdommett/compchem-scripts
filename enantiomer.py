#!/usr/bin/env python
# simple script which reflects a molecule in the z axis
from sys import argv

mol_file=open(argv[1]).read().splitlines()
outfile=open(argv[1][:-4]+"_enantiomer.xyz","w")

natoms=int(mol_file[0].split()[0])
outfile.write("{}\n\n".format(natoms))
for line in mol_file[2:]:
    symb,x,y,z=line.split()
    outfile.write("{0:<2} {1:>13.9s} {2:>13.9s} {3:>13.9f}\n".format(symb,x,y,-float(z)))
outfile.close()
