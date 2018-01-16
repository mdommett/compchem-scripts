#!/usr/bin/env python2
import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element

infile = open(argv[1],"r").read().splitlines()

for line in infile:
	if "NAtoms=" in line:
		natoms = int(line.split()[3])
		break

def get_charges(inputfile,natoms):

	q = []
	count= []
	for linenumber,line in enumerate(infile):
		if line == " ESP charges:":
			count.append(int(linenumber))
	for line in count:

		charges = [i.split() for i in infile[(line+2):(line+2+natoms)]]
		q.append(charges)

	return q

ESP_charges = (get_charges(infile,natoms))

for set in range(len(ESP_charges)):
	for line in ESP_charges[set]:
		print "{0:<2}  {1:>3}  {2:>8.5f}".format(int(line[0]), str(line[1]), float(line[2]))
