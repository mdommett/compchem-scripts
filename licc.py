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
import sys
import copy

def interpolate(initial,final,diff,steps):
    for step in range(steps+1):
        new_coords=copy.deepcopy(initial)
        for atom in range(len(initial)):
            new_coords[atom].x=initial[atom].x+(diff[atom,0]/steps)*step
            new_coords[atom].y=initial[atom].y+(diff[atom,1]/steps)*step
            new_coords[atom].z=initial[atom].z+(diff[atom,2]/steps)*step
        ef.write_xyz(str('step_{}.xyz'.format(step)),new_coords)

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--initial", help="Initial .xyz file", type=str)
    parser.add_argument("-f","--final", help="Final.xyz file", type=str)
    parser.add_argument("-n","--nsteps", help="Number of steps", type=int)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    init_xyz=rf.read_pos(args.initial)
    final_xyz=rf.read_pos(args.final)

    init=np.array([(i.x,i.y,i.z) for i in init_xyz])
    final=np.array([(i.x,i.y,i.z) for i in final_xyz])
    diff=final-init

    interpolate(init_xyz,final_xyz,diff,args.nsteps)
