#!/usr/bin/env python
"""
Calculate the adiabatic reorganisation energy from two save_geometries

usage: $adiabatic_reorganisation_energy.py ground_state_minimum.out excited_state_mininum.out
"""
import sys
gs_string="SCF Done:"
es_string="Total Energy"
def scf_energy(infile):
    gs=[]
    for line in infile:
        if gs_string in line:
            gs.append(float(line.split()[4]))
    return gs[-1]

def es_energy(infile):
    es=[]
    for line in infile:
        if es_string in line:
            es.append(float(line.split()[4]))
    return es[-1]

gs=open(sys.argv[1],"r").read().splitlines()
es=open(sys.argv[2],"r").read().splitlines()
gs_s0=scf_energy(gs)
gs_s1=es_energy(gs)
es_s0=scf_energy(es)
es_s1=es_energy(es)

reorg_energy_au=(gs_s1-es_s1)+(es_s0-gs_s0)
reorg_energy_ev=reorg_energy_au*27.2114
print("Lambda = {0:.2f} eV".format(reorg_energy_ev))
