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

def dihedral_angle(p0,p1,p2,p3):
    # q0 = p1-p0 vector connecting p0 to p1
    q0 = p1 - p0
    q1 = p2 - p1
    q2 = p3 - p2

    q0q1=np.cross(q0,q1) # vector normal to plane p0p1p2
    q1q2=np.cross(q1,q2) # normal to plane p1p2p3

    n0=q0q1/magnitude(q0q1) # unit normal
    n1=q1q2/magnitude(q1q2)

    u0=n1 # orthogonal unit vectors
    u2=q1/magnitude(q1)
    u1=np.cross(u0,u2)

    cos_theta=np.dot(n0,u0)
    sin_theta=np.dot(n0,u1)
    theta=-np.arctan2(sin_theta,cos_theta)
    return np.degrees(theta)

def magnitude(vector):
    return np.linalg.norm(vector)

if __name__=='__main__':
    geom_file=rf.read_pos(argv[1])
    geom=coordinate_matrix(geom_file)
    A=geom[int(argv[2])-1]
    B=geom[int(argv[3])-1]
    C=geom[int(argv[4])-1]
    D=geom[int(argv[5])-1]

    dihedral_angle=dihedral_angle(A,B,C,D)
    print(dihedral_angle)
