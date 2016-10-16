#!/usr/bin/env python

from sys import argv,exit
from collections import defaultdict
import numpy as np

""" This script extracts atoms from a template to create 2 files"""

template = open(argv[1],"r").read().splitlines()
qm_atoms = open(argv[2],"r").read().splitlines()
outfile_qm = open(argv[3],"w")
outfile_mm = open(argv[4],"w")

for atom in qm_atoms:
    atom_no = int(atom)+1
    outfile_qm.write("{0}\n".format(template[atom_no]))

int_qm_atoms = []
for atom in qm_atoms:
    atm_no = int(atom)
    int_qm_atoms.append(atm_no)

count = 1
for line in template[2:]:

    if count not in int_qm_atoms:
        outfile_mm.write("{0}\n".format(line))
        count +=1
    else:
        count +=1

