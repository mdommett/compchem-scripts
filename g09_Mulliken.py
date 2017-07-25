#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element

infile = open(argv[1],"r").read().splitlines()

for line in infile:
	if "NAtoms=" in line:
		natoms = int(line.split()[1])
		break

def get_charges(inputfile,natoms):

    q = []

    for i,j in enumerate(infile):
        if j == " Mulliken charges:":
            count = int(i)
            break
    charges = infile[(count+2):(count+2+natoms)]
    summed = infile[count+natoms+2].split()
    for line in charges:
        line_charge= line.split()
        q.append(float(line_charge[2]))
    sumd=0
    for i in q:
        sumd +=float(i)
    return q

charges = (get_charges(infile,natoms))

for i in charges:
	print "{:>9.6f}".format(i)
