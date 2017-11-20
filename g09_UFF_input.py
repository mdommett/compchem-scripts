#!/usr/bin/env python

import numpy as np
import read_file as rf
import edit_file as ef
import handle_atoms as ha
import dimer_select as ds
from math import sqrt
import argparse
import sys

def make_clust(atoms, clust_rad, max_bl, refx, refy, refz,natoms):
    # atoms within the sphere of rad clust_rad
    seed_atoms = []
    for atom in atoms:
        if atom.dist(refx,refy,refz) < clust_rad:
            seed_atoms.append(atom)
    ef.write_xyz("seedatoms.xyz",seed_atoms)

    # atoms in the cluster (seed_atoms + atoms to complete molecules)
    clust_atoms = []
    for i,atom in enumerate(seed_atoms):
        if atom not in clust_atoms:
            mol2Add = ha.select(max_bl, atoms, atoms.index(atom))
            if len(mol2Add)==natoms:
                for atom2Add in mol2Add:
                    clust_atoms.append(atom2Add)
    return clust_atoms


parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input .xyz file",type=str)
parser.add_argument("-b", "--bond", help="Maximum length (in unites of input file) that qualifies as a bond",
                        default=1.6, type=float)
parser.add_argument("-r", "--radius", help="Radius of cluster",type=float,default=10)
user_input = sys.argv[1:]
args = parser.parse_args(user_input)
cluster=rf.read_xyz(args.input)[-1]
cluster_centroid=ha.find_centroid(cluster)
distances=[ds.vector_distance((cluster_centroid[0],cluster_centroid[1],cluster_centroid[2],atom.x,atom.y,atom.z)) for atom in cluster]
central_atom=cluster[distances.index(min(distances))]
central_molecule=ha.select(args.bond,cluster,distances.index(min(distances)))
natoms=len(central_molecule)
for i in central_molecule:
    cluster.remove(i)
    cluster.insert(0,i)
c_m_centroid=ha.find_centroid(central_molecule)
final_cluster=make_clust(cluster,args.radius,args.bond,c_m_centroid[0],c_m_centroid[1],c_m_centroid[2],natoms)
print c_m_centroid
print ha.find_centroid(final_cluster)
ef.write_xyz(args.input[:-4]+"-cluster.xyz",final_cluster)
######################################################################

jobname=args.input[:-4]+"_UFF"
comfile = open(jobname+".com","w")

chk = "%chk={0}.chk".format(jobname)
nproc = "%nproc=4"
mem = "%mem=28GB"
input_line = "#p ONIOM(wb97xd/6-31G*:UFF=QEq)=(EmbedCharge) force nosymm"
comfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n".format(chk,nproc,mem,input_line))
for atom in final_cluster[:natoms]:
    atomStr = "{:>6} 0 {:10.6f} {:10.6f} {:10.6f} H \n".format(atom.elem, atom.x, atom.y, atom.z)
    comfile.write(atomStr)
for atom in final_cluster[natoms:]:
    atomStr = "{:>6} -1 {:10.6f} {:10.6f} {:10.6f} M \n".format(atom.elem, atom.x, atom.y, atom.z)
    comfile.write(atomStr)
#comfile.write("\n--link1--\n")
#input_line="#p ONIOM(wb97xd/6-311++G** td=(nstates=3):UFF=QEq)=(EmbedCharge) geom=check guess=read nosymm"
#comfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n\n".format(chk,nproc,mem,input_line))
