#!/usr/bin/env python

import sys
import read_file as rf

atoms = rf.read_xyz(sys.argv[1])[-1]
outfile = open(sys.argv[1][:-4]+".com","w")

chk = "%chk={0}.chk".format(sys.argv[1][:-4])
nproc = "%nproc=4"
mem = "%mem=28GB"
input_line = "#p wb97xd 6-311++G** td=(nstates=4) nosymm"
###############
outfile.write("{0}\n{1}\n{2}\n{3}\n\n Title \n\n0 1\n".format(chk,nproc,mem,input_line))
for atom in atoms:
    atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
    outfile.write(atomStr)
outfile.write("\n")
