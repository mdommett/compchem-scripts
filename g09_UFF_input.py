#!/usr/bin/env python

import numpy as np
import read_file as rf
import edit_file as ef
import handle_atoms as ha
import dimer_select as ds
from math import sqrt
import argparse
import sys
from operator import itemgetter
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input .xyz file",type=str)
parser.add_argument("-b", "--bond", help="Maximum length (in unites of input file) that qualifies as a bond",
                        default=1.6, type=float)
parser.add_argument("-r", "--radius", help="Radius of cluster",type=float)
user_input = sys.argv[1:]
args = parser.parse_args(user_input)

cell=rf.read_xyz(args.input)[-1]
cluster=ha.make_cluster(cell,args.radius,args.bond)
cluster_centroid=ha.find_centroid(cluster)
molecules=ds.make_molecules(cluster,args.bond)
natoms=len(molecules[0])
unordered=[]
for molecule in molecules:
    centroid=ha.find_centroid(molecule)
    new_list=[molecule,ds.vector_distance((centroid[0],centroid[1],centroid[2],cluster_centroid[0],cluster_centroid[1],cluster_centroid[2]))]
    unordered.append(new_list)
ordered=sorted(unordered,key=itemgetter(1))
final_molecules=[i[0] for i in ordered]
final_cluster=[item for sublist in final_molecules for item in sublist]
ef.write_xyz(args.input[:-4]+"-cluster.xyz",final_cluster)
jobname=args.input[:-4]+"_UFF
comfile = open(jobname+".com","w")

chk = "%chk={0}.chk".format(jobname)
nproc = "%nproc=4"
mem = "%mem=28GB"
input_line = "#p ONIOM(wb97xd/6-31Gc*UFF=QEq)=(EmbedCharge) nosymm"
comfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n".format(chk,nproc,mem,input_line))
for atom in final_cluster[:natoms]:
    atomStr = "{:>6} 0 {:10.6f} {:10.6f} {:10.6f} H \n".format(atom.elem, atom.x, atom.y, atom.z)
    comfile.write(atomStr)
for atom in final_cluster[natoms:]:
    atomStr = "{:>6} -1 {:10.6f} {:10.6f} {:10.6f} M \n".format(atom.elem, atom.x, atom.y, atom.z)
    comfile.write(atomStr)
comfile.write("\n--link1--\n")
input_line="#p ONIOM(wb97xd/6-311++G** td=(nstates=3):UFF=QEq)=(EmbedCharge) geom=check guess=read nosymm"
comfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n\n".format(chk,nproc,mem,input_line))


ef.write_xyz("test.xyz",cluster)
