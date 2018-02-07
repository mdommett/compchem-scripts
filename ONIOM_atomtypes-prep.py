#!/usr/bin/env python

from sys import argv,exit
atoms=open("atoms","r").read().splitlines()
types=open("atomtypes","r").read().splitlines()

with open("atomtypes-complete","w") as atomtypes_tmp:
    for i in range(len(atoms)):
        atomtypes_tmp.write("{}-{}\n".format(atoms[i],str(types[i].split()[0])))
    atomtypes_tmp.close()
