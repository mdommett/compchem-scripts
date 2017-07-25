#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys


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

def test_range(chrom1,chrom2):
    print "Unit 1: {0} atoms".format(len(np.unique(chrom1)))
    print "Unit 2: {0} atoms".format(len(np.unique(chrom2)))
    if len(np.unique(chrom1+chrom2)) != natoms:
        print "You have made a mistake in the molecule input!\nThe number of atoms specified is not equal to the total.\nPlease fix the molecule input and retry."
        missing = sorted(set(range(1, natoms)).difference(np.unique(chrom1+chrom2)))
        print "The missing atom number(s) is: {0}".format(missing)
        sys.exit()
    if len(chrom1)+len(chrom2) != natoms:
        print "The total number of atoms is {0}, but there are {1} atoms in the two units!\nPlease fix the molecule input and retry.".format(natoms,len(chrom1)+len(chrom2))
        sys.exit()


molecule = open("molecule","r")
molecule_lines = molecule.read().splitlines()
natoms = int(molecule_lines[0])
chrom1 = rangeexpand(molecule_lines[1])
chrom2 = rangeexpand(molecule_lines[2])
test_range(chrom1,chrom2)

state1 = open(argv[1],"r").read().splitlines()
state1_charges = [float(i.split()[2]) for i in state1]
state2= open(argv[2],"r").read().splitlines()
state2_charges = [float(i.split()[2]) for i in state2]

total_chrom1 = 0
total_chrom2 = 0   

for i in chrom1:
    total_chrom1 += (state2_charges[i-1]-state1_charges[i-1])
for i in chrom2:
    total_chrom2 += (state2_charges[i-1]-state1_charges[i-1])
print "Unit 1: {}\nUnit 2: {}\n".format(total_chrom1,total_chrom2)
