#!/usr/bin/env python

from sys import argv,exit

def parse_energies(infile):
    E,f=[],[]
    for line in infile:
        if " Excited State" in line:
            E.append(float(line.split()[4]))
            f.append(float(line.split()[8][2:]))
    if not E or not f:
        exit("Could not find any excitation energies")
    return E[-1],f[-1]

def ev_to_cm(value):
    return value/8065.54429
def radiative_rate(E,f):
    return (f*(ev_to_cm(E)**2))/(1.499)

if __name__=='__main__':
    infile=open(argv[1],"r").read().splitlines()
    E,f=parse_energies(infile)
    K=radiative_rate(E,f)
    print("E = {:.2f} eV".format(E))
    print("f = {:.3f}   ".format(f))
    print("Calculating the radiative decay rate...")
    print("-----------------")
    print("K = {:.3g} s^-1".format(K))
    print("1/K = {:.3g} s".format(1/K))
