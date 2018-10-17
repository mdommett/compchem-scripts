#!/usr/bin/env python
import numpy as np
from sys import argv,exit
from sklearn.metrics.pairwise import pairwise_distances

## Reports the tau paramter (https://en.wikipedia.org/wiki/Geometry_index)
# usage:
# tau_paramter file.xyz central_atom_number
# Returns:
# tau parameter
def file_to_matrix(infile):
    xyz=np.zeros((len(infile),3))
    for i in range(len(infile)):
        xyz[i,0]=infile[i].split()[1]
        xyz[i,1]=infile[i].split()[2]
        xyz[i,2]=infile[i].split()[3]
    return xyz

def angle(point1,central,point2):
    ba=point1-central
    bc=point2-central
    cosine_angle=np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


infile=open(argv[1],'r').readlines()

try:
    natoms=int(infile[0].split()[0])
    if natoms!=len(infile[2:]):
        exit("Number of atoms doesn't match number specified in xyz")
except ValueError:
    "Cannot read number of atoms!"
xyz=file_to_matrix(infile[2:])
central_atom=int(argv[2])-1
distance_matrix=pairwise_distances(xyz)
coordinating=np.argsort(distance_matrix[central_atom,:])[1:6]
print("Coordinating Atom Numbers: {}".format(coordinating+1))
angles=[]
for i,p1 in enumerate(coordinating):
    for j,p2 in enumerate(coordinating[i:]):
        if p1!=p2:
            angles.append(angle(xyz[p1],xyz[central_atom],xyz[p2]))
sorted_angles=sorted(angles)
beta=sorted_angles[-1]
alpha=sorted_angles[-2]
tau=(beta-alpha)/60
print(tau)
