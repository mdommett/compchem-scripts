#!/usr/bin/env python3
"""
script to parse sh.log file from SOC dynamics and extract the SOC values. Usage:
python3 parse_sh_SOC.py nsinglets ntriplets singlet_coupling_state

nsinglets and ntriplets are the total number of singlet and triplet states (MCH rep.) used in the dynamics. nsinglets must include S0

singlet_coupling_state is the singlet state you want the coupling to:
s0=0
s1=1
s2=2
etc

eg parse_sh_SOC.py 4 4 1 >couplings

triplet components are added up

Michael Dommett
26/6/2018
"""

import numpy as np
from sys import argv,exit

nsinglets=int(argv[2])
ntriplets=int(argv[3])
singlet_coupling_state=int(argv[4])
nstates=nsinglets+3*ntriplets
couplings=[]
with open(argv[1]) as infile:
    for line in infile:
        if "#DEBUG# Couplings (cm-1) =" in line:
            coupling_matrix=np.zeros((nstates,nstates))
            for i in range(nstates):
                coupling_line=next(infile).split()
                for j,x in enumerate(coupling_line):
                    a, b = x.strip('()').split(',')
                    coupling_matrix[i,j]=a

            couplings.append(coupling_matrix)
a=[nsinglets,nsinglets+ntriplets,nsinglets+ntriplets+ntriplets]
states=[]
states.append(a)
for i in range(1,ntriplets):
    states.append([j+i for j in a])
for matrix in couplings:
    printlist=[np.linalg.norm(matrix[singlet_coupling_state,state]) for state in states]
    print(*printlist,sep = "    ")
