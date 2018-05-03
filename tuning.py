#!/usr/bin/env python
"""
Calculates the J^2 value for a set of N,N+1 and N-1 electron systems

J^2(y) = (Ey_homo(N) + IPy(N))**2 + (Ey_homo(N+1) + IPy(N+1))**2

usage:
tuning.py N.out N+1.out N-1.out

Michael Dommett
May 2018
"""
from sys import argv

def parse_energy(infile):
    E=[]
    for line in infile:
        if "SCF Done" in line:
            E.append(line.split()[4])
    return float(E[-1])

def parse_orbital(infile):
    orbital_energies=[]
    for line in infile:
        if "Alpha  occ. eigenvalues" in line:
            #we dont know how many orbital energies there are,
            # so we just take a list
            orbital_energies.append((line.split()[4:]))
    # orbital energies is a list of lists, so we want the final
    # value of the final list
    return float(orbital_energies[-1][-1])

def calculate_IP(neutral,cation):
    return parse_energy(cation)-parse_energy(neutral)

# N electron system
N=open(argv[1]).read().splitlines()
# N+1 electron system (anion)
N_p_1=open(argv[2]).read().splitlines()
# N-1 electron system (cation)
N_m_1=open(argv[3]).read().splitlines()

Ehomo_N=parse_orbital(N)
IP_N=calculate_IP(N,N_m_1)
Ehomo_N_p_1=(parse_orbital(N_p_1))
IP_N_p_1=calculate_IP(N_p_1,N)

J=(Ehomo_N+IP_N)**2+(Ehomo_N_p_1+IP_N_p_1)**2

print("{:.4g}".format(J))
