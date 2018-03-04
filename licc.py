#!/usr/bin/env python
"""
Linear interpolation of cartesian coordinates

Usage:
python -i init_geom -f final_geom -n no. of steps

output: licc_*.xyz
"""
import read_file as rf
import edit_file as ef
import numpy as np
import argparse

def interpolate(initial,final,steps):
    for step in range(len(steps)+1):
        new_coords==init_xyz
        for atom in range(len(init)):
            new_coords[atom].x=init_xyz[atom].x-(diff[atom,0]/steps)*step
            new_coords[atom].y=init_xyz[atom].y-(diff[atom,1]/steps)*step
            new_coords[atom].z=init_xyz[atom].z-(diff[atom,2]/steps)*step
    ef.write_xyz('_{}.xyz'.format(step),new_coords)

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--initial" help="Initial .xyz file", type=str)
    parser.add_argument("-f","--final" help="Final.xyz file", type=str)
    parser.add_argument("-n","--nsteps" help="Number of steps", type=int)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    interpolate(args.initial,args.final,args.nsteps)
