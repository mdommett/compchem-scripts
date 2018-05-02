#!/usr/bin/env python
"""
Calculates the distance between centroids of two molecules in one xyz file

THE MOLECUES MUST BE IN THE .xyz FILE ONE AFTER ANOTHER
e.g for a water dimer:
----------
6

H x1 y1 z1
H x2 y2 z2
O x3 y3 z3
H x4 y4 z6
H x5 y5 z6
O x6 y6 z6
----------

usage: dimer_centroids.py dimer.xyz
output: distance in units of xyz file

Michael Dommett
May 2018
"""
import read_file as rf
import dimer_select as ds
import handle_atoms as ha
from sys import argv

def centroid_distance(dimer):
    """Calculates the centroid distances between two molecules of a dimer"""
    cent_1=ha.find_centroid(dimer[0:int(len(dimer)/2)])
    cent_2=ha.find_centroid(dimer[int(len(dimer)/2):len(dimer)])
    return ds.vector_distance(cent_1[0],cent_1[1],cent_1[2],cent_2[0],cent_2[1],cent_2[2])

dimer=rf.read_pos(argv[1])
print("{:.3f}".format(centroid_distance(dimer)))
