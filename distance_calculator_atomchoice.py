#!/usr/bin/env python

import numpy as np
from sys import argv,exit
from itertools import chain
from math import sqrt
from collections import defaultdict
#import matplotlib.pyplot as plt
import periodic as pd

####
# Usage: $distance_calculator.py $1.xyz molecule
#
#
# The molecule file is structured:
#
# number of atoms to explore
# reference atom numbers in xyz file
# atom numbers in xyz to get the distance to
# specifc atom symbols you want to get the distances for, space separated
#
# Example:
#
# 1037
# 9-34, 35, 28, 21, 63
# 100-1032, 1036, 1-8
# C H N O
# The output is returned to the terminal, in descending
# order of distances:
# (ref atom no.,atom no.): distance (in units of .xyz file)
#
# Enjoy!
#
# Michael Dommett, March 2017
# m.dommett@qmul.ac.uk
#####

#### Important Functions ####

# Split up the list of numbers

def group_to_range(group):
    group = ''.join(group.split())
    sign, g = ('-', group[1:]) if group.startswith('-') else ('', group)
    r = g.split('-', 1)
    r[0] = sign + r[0]
    r = sorted(int(__) for __ in r)
    return range(r[0], 1 + r[-1])

# Expand and sort the list of numbers
def rangeexpand(txt):
    ranges = chain.from_iterable(group_to_range(__) for __ in txt.split(','))
    return sorted(set(ranges))

# Extract x,y,z positions of given line

def get_coordinates(coordinates):

    x = float(coordinates.split()[1])
    y = float(coordinates.split()[2])
    z = float(coordinates.split()[3])
    return x,y,z

# Calculate distance between 3 (x,y,z) points

def vector_distance((x1,y1,z1,x2,y2,z2)):
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2-z1)**2)
    return dist

def bin_width(dataset):

    q75,q25 = np.percentile(dataset,[75,25])
    IQR = q75-q25

    return 2*(IQR/(len(dataset)**(1./3)))

def number_of_bins(dataset):

    return round((max(dataset)-min(dataset))/bin_width(dataset))

def plot_histogram(dataset):
    weights = np.ones_like(dataset)/len(dataset)
    n, bins, patches = plt.hist(
    dataset,
    #normed=1,
    alpha=0.6,weights=weights)
    return

### Sript ###

# Load input parameters

xyz = open(argv[1],'r').readlines()
molecule = open("molecule",'r').readlines()
natoms= int(xyz[0].split()[0])

# Set input definitions

ref_atoms = rangeexpand(molecule[1])
cluster_atoms = rangeexpand(molecule[2])
atoms_to_search = [i for i in molecule[3].split()]

### Check the atomic inputs match eachother

if natoms != int(molecule[0].split()[0]):
    print "Molecule number of atoms not equal to xyz number of atoms"
    exit
if natoms != len(xyz[2:]):
    print "xyz file has {0} atom coordinates but {1} atoms specified in preamble".format(len(xyz[2:]),natoms)


## Calculate distances

# Initiate dictionary
distances = defaultdict(list)

# Populate dictionary

for i in ref_atoms:
    for k in [int(l) for l in cluster_atoms if xyz[1+l].split()[0] in atoms_to_search]:
        distances[(i,k,xyz[1+k].split()[0])] =(vector_distance(get_coordinates(xyz[1+i])+get_coordinates(xyz[1+k])))

# Print dictionary sorted by distances

d_view = [ (v,k) for k,v in distances.iteritems() ]
d_view.sort(reverse=True)
for v,k in d_view:
    print "{0} - {1} ({2}) : {3:.3f}".format(k[0],k[1],k[2],v)
#plot_histogram([v for k,v in distances.iteritems()])
#plt.show()

