#!/usr/bin/env python

from sys import argv

inputfile = open(argv[1],"r").read().splitlines()
outputfile = open(argv[2],"w")
natoms = int(inputfile[0])


frozen_atoms = raw_input("Enter atom numbers to be frozen: ").split()
frozen_atoms_list = [int(a)+1 for a in frozen_atoms]
count = 2
outputfile.write("{0}\n\n".format(natoms))

for atom_line in inputfile[2:natoms+2]:
    atom_line_split = atom_line.split()
    if count in frozen_atoms_list:
        atom_line_split.insert(1,-1)
    else:
        atom_line_split.insert(1,0)
    for i in atom_line_split:
        outputfile.write("{0}   ".format(i))
    outputfile.write("\n")
    count += 1