#!/usr/bin/env python
import read_file as rf
from sys import argv
import numpy as np
import edit_file as ef

"""
ref: https://pubs.acs.org/doi/suppl/10.1021/acs.jpca.5b06639/suppl_file/jp5b06639_si_001.pdf
    A
    |
 C--B--D
"""

def coordinate_matrix(atoms):
    coords=np.zeros((len(atoms),3))
    for i in range(len(atoms)):
        coords[i,0]=atoms[i].x
        coords[i,1]=atoms[i].y
        coords[i,2]=atoms[i].z
    return coords

def magnitude(vector):
    return np.linalg.norm(vector)

if __name__=='__main__':
    geom_file=rf.read_pos(argv[1])
    geom=coordinate_matrix(geom_file)
    A=int(argv[2])-1
    B=int(argv[3])-1
    C=int(argv[4])-1
    D=int(argv[5])-1

    BA=geom[A]-geom[B]
    BC=geom[C]-geom[B]
    BD=geom[D]-geom[B]

    angle=np.arccos(np.dot(np.cross(BC,BD),BA)/np.dot(magnitude(np.cross(BC,BD)),magnitude(BA)))
    pyr_angle_radians=np.pi/2 - angle
    pyr_angle=np.degrees(pyr_angle_radians)
    print(pyr_angle)
