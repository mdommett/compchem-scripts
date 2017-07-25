#!/usr/bin/env python

"""----------------
This script diabatizes two non-adiabatic states to produce a diabatic Hamiltonion,
the off-diagonal elements of which are the couplings J between the two electronic
states.

The diabatization scheme used herein is proposed by Troisi et. al PRL 114, 026402 (2015)
 (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.026402).
 The major detail can be found in the supplementary information,
 where they argue that the exciton coupling J can be
 computed by diabatization of the adiabatic s1 and s2 states of a dimer.
 The best 'diabatic' transformation is yielded by the matrix which minimizes
 the difference between the transition dipole moments of the s1 and s2 states
 of the dimer and the s1 states of the constituent monomers.

 Requirements: TDM and energies of the s1 and s2 states of the dimer,
                TDM of the s1 states for the constinuent monomers of the dimer

Usage:

./diabatization.py dimer_file.out monomer-1.out monomer-2.out

The three G09 output files must be calculated with the "nosymm" option

For the dimer, the first two excited states must be calculated
For the two constituent monomers, in the dimer config, the first excited states
must be calculated

Michael Dommett
June 2017
m.dommett@qmul.ac.uk
----------------
"""
import numpy as np
from sys import argv
import pandas as pd

def get_tdm(file,state):
    for linenumber,line in enumerate(file):
        if line==" Ground to excited state transition electric dipole moments (Au):":
            stateline=1+int(state)+linenumber
            break
    return np.matrix(file[stateline].split()[1:4],dtype=float)

def get_energy(file,state):
    for line in file:
        if " Excited State   {}".format(state) in line:
            energy=float(line.split()[4])
            break
    return energy


def diabatize(TDMs1,TDMs2,tdmA,tdmB,Es1,Es2):
    TDMdim=np.concatenate((TDMs1,TDMs2))
    TDMmon=np.concatenate((tdmA,tdmB))
    M=np.dot(TDMdim,TDMmon.transpose())

    U,s,Vt= np.linalg.svd(M)

    C=(np.dot(U,Vt)).transpose()

    E=np.matrix(([Es1,0],[0,Es2]))

    H=np.dot(np.dot(C,E),C.transpose())

    return H


dimer=open(argv[1],"r").read().splitlines()
monA=open(argv[2],"r").read().splitlines()
monB=open(argv[3],"r").read().splitlines()

# Dimer data (s1 and s2)
TDMs1=get_tdm(dimer,1)
TDMs2=get_tdm(dimer,2)
Es1=get_energy(dimer,1)
Es2=get_energy(dimer,2)

# Monomer data (s1)
tdmA=get_tdm(monA,1)
tdmB=get_tdm(monB,1)
H=diabatize(TDMs1,TDMs2,tdmA,tdmB,Es1,Es2)
J=H[0,1].round(3)
print "Diabatic Hamiltonian H:\n{}\n".format(H)
print "J = {} eV".format(J)
