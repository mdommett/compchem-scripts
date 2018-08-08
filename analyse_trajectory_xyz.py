#!/usr/bin/env python3

import numpy as np
from sys import argv,exit
import argparse
import read_file as rf

def magnitude(vector):
    return np.linalg.norm(vector)
def pyramidal_angle(geom,A,B,C,D):
    A=int(A)-1
    B=int(B)-1
    C=int(C)-1
    D=int(D)-1

    BA=geom[A]-geom[B]
    BC=geom[C]-geom[B]
    BD=geom[D]-geom[B]

    angle=np.arccos(np.dot(np.cross(BC,BD),BA)/np.dot(magnitude(np.cross(BC,BD)),magnitude(BA)))
    pyr_angle_radians=np.pi/2 - angle
    pyr_angle=np.degrees(pyr_angle_radians)

    return pyr_angle

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dihedral", help="Dihedral angle between atoms A,B,C,D", nargs=4,type=list)
    parser.add_argument("-a","--angle",help="Angle between atoms A,B,C",nargs=3,type=list)
    parser.add_argument("-p","--pyramidal",help="pyramidal angle between bond AB plane CD",
    nargs=4,type=list)
    parser.add_argument("input", help="Input files",type=str,nargs='*')
    args = parser.parse_args(user_input)
    
    geoms=rf.read_xyz(input)
