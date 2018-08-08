#!/usr/bin/env python3

from sys import argv
import numpy as np
from cryspy.io import read_file as rf

def coordinate_matrix(atoms):
    coords=np.zeros((len(atoms),3))
    for i in range(len(atoms)):
        coords[i,0]=atoms[i].x
        coords[i,1]=atoms[i].y
        coords[i,2]=atoms[i].z
    return coords

def bond_length(A,B):
    return np.linalg.norm(B-A)

if __name__=='__main__':
    geom_file=rf.read_pos(argv[1])
    geom=coordinate_matrix(geom_file)
    A=geom[int(argv[2])-1]
    B=geom[int(argv[3])-1]
    bond_length=bond_length(A,B)
    print("{:.3g}".format(bond_length))
