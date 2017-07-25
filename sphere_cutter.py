#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element
from math import sqrt

def group_to_range(group):
    group = ''.join(group.split())
    sign, g = ('-', group[1:]) if group.startswith('-') else ('', group)
    r = g.split('-', 1)
    r[0] = sign + r[0]
    r = sorted(int(__) for __ in r)
    return range(r[0], 1 + r[-1])

# Expand and sort the list of numbers
def rangeexpand(txt):
    ranges = chain.from_iterable(group_to_range(__) for __ in txt.split(','))
    return sorted(set(ranges))

def get_coordinates(coordinates):
    symbol=float(coordinates.split()[0])
    x = float(coordinates.split()[1])
    y = float(coordinates.split()[2])
    z = float(coordinates.split()[3])
    return x,y,z

def vector_distance((x1,y1,z1,x2,y2,z2)):
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2-z1)**2)
    return dist

def centre_of_mass(listofatoms):
    Xnumerator,Ynumerator,Znumerator,Xdenominator,Ydenominator,Zdenominator=0,0,0,0,0,0
    for index in listofatoms:
        line = xyz[index+1]
        symbol,x,y,z = str(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])

        Xnumerator += element(symbol).mass*x
        Ynumerator += element(symbol).mass*y
        Znumerator += element(symbol).mass*z
        Xdenominator += element(symbol).mass
        Ydenominator += element(symbol).mass
        Zdenominator += element(symbol).mass

    Xcm = Xnumerator/Xdenominator
    Ycm = Ynumerator/Ydenominator
    Zcm = Znumerator/Zdenominator

    return Xcm,Ycm,Zcm

# Load input parameters

xyz = open(argv[1],'r').readlines()
molecule = open("molecule",'r').readlines()
natoms= int(xyz[0].split()[0])
radius=float(raw_input("Radius of sphere in Angstroms: "))
outfile_answer = raw_input("Do you want to write an xyz file for the output? [Y/N] ")
# Set input definitions

ref_atoms = rangeexpand(molecule[1])
cluster_atoms = rangeexpand(molecule[2])
#atoms_to_search = [i for i in molecule[3].split() if molecule[3]]

if natoms != int(molecule[0].split()[0]):
    print "Molecule number of atoms not equal to xyz number of atoms"
    exit
if natoms != len(xyz[2:]):
    print "xyz file has {0} atom coordinates but {1} atoms specified in preamble".format(len(xyz[2:]),natoms)

com=centre_of_mass(ref_atoms)
sphere_atoms=[]
for j in cluster_atoms:
    x2,y2,z2=xyz[1+j].split()[1:]
    if vector_distance((com[0],com[1],com[2],float(x2),float(y2),float(z2))) <= radius:
        sphere_atoms.append(j)

if outfile_answer == "Y":
    outfile = open(raw_input("Name of outfile: "),"w")
    outfile.write("{}\n\n".format(len(ref_atoms)+len(sphere_atoms)))
    for i in ref_atoms:
        outfile.write("{}".format(xyz[1+i]))
    for i in sphere_atoms:
        outfile.write("{}".format(xyz[1+i]))

comfile=open(argv[2],'r').readlines()
com_new=open(argv[2]+"_new","w")

for i in range(len(comfile)):
    if int(i)-8 not in sphere_atoms:
        com_new.write("{}".format(comfile[i-1]))
    else:
        line = comfile[i-1].split()
        com_new.write("{0:<2}  0 {1:>13.9f} {2:>13.9f} {3:>13.9f} {4}\n".format(line[0],float(line[2]),float(line[3]),float(line[4]),line[5]))
