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

def test_range(metal,chrom1,chrom2):
    print "Unit 1: {0} atoms".format(len(np.unique(metal)))
    print "Unit 2: {0} atoms".format(len(np.unique(chrom1)))
    print "Unit 3: {0} atoms".format(len(np.unique(chrom2)))
    print "Total: {} atoms".format(len(np.unique(metal+chrom1+chrom2)))
    if len(np.unique(chrom1+chrom2+metal)) != natoms:
        print "You have made a mistake in the molecule input!\nThe number of atoms specified is not equal to the total.\nPlease fix the molecule input and retry."
        missing = sorted(set(range(1, natoms+1)).difference(np.unique(chrom1+chrom2+metal)))
        print "The missing atom number(s) is: {0}".format(missing)
        sys.exit()
    if len(chrom1)+len(chrom2)+len(metal) != natoms:
        print "The total number of atoms is {0}, but there are {1} atoms in the two units!\nPlease fix the molecule input and retry.".format(natoms,len(chrom1)+len(chrom2))
        sys.exit()


molecule = open(argv[2],"r")
molecule_lines = molecule.read().splitlines()
natoms = int(molecule_lines[0])
metal = rangeexpand(molecule_lines[1])
rotaxane = rangeexpand(molecule_lines[2])
skeleton = [i for i in range(1,natoms+1) if i not in metal+rotaxane ]
print range(1,natoms)

test_range(metal,rotaxane,skeleton)

NBOfile = open(argv[1],"r").read().splitlines()
NBO_charges = [float(i.split()[2]) for i in NBOfile]

total_metal=0
total_rotaxane = 0
total_skeleton= 0

for i in metal:
    total_metal += (NBO_charges[i-1])
for i in rotaxane:
    total_rotaxane += (NBO_charges[i-1])
for i in skeleton:
    total_skeleton += (NBO_charges[i-1])
print "Metal: {}\nRotaxane: {}\nSkeleton:{}".format(total_metal,total_rotaxane,total_skeleton)
print "Total: {}".format(total_skeleton+total_metal+total_rotaxane)
