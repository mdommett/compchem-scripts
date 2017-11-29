#!/usr/bin/env python

from sys import argv,exit
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-at","--AtomTypes",help="atom types file", type=str)
parser.add_argument("-ac","--AtomicCharges",help="atom charges file", type=str)
parser.add_argument("-tmp","--TmpAtoms",help="active atoms xyz ",type=str)
parser.add_argument("-f","--FrozenAtoms",help="frozen atoms xyz",type=str)
user_input = argv[1:]
args=parser.parse_args(user_input)

AtomTypes=open(args.AtomTypes,"r").read().splitlines()
AtomicCharges=open(args.AtomicCharges,"r").read().splitlines()
TmpAtoms=open(args.TmpAtoms,"r").read().splitlines()[2:]
FrozenAtoms=open(args.FrozenAtoms,"r").read().splitlines()[2:]

if len(AtomTypes)!=len(AtomicCharges) or len(TmpAtoms)+len(FrozenAtoms)!=len(AtomTypes):
	exit("Input is incorrect - there are {} atom types, {} charges, {} tmp atoms \
	and {} frozen atoms.\nFix and try again").format(len(AtomTypes), \
	len(AtomicCharges),len(TmpAtoms),len(FrozenAtoms))
else:
	with open("tmp.com","w") as tmpcom:
		tmpcom.write("{}\n\n".format(len(TmpAtoms)))
		for i in range(len(TmpAtoms)):
			tmpcom.write("{}\n".format(TmpAtoms[i]))
		tmpcom.close()
	with open("frozen.xyz","w") as frznxyz:
		frznxyz.write("{}\n\n".format(len(FrozenAtoms)))
		for i in range(len(FrozenAtoms)):
			frznxyz.write("{}\n".format(FrozenAtoms[i]))
		frznxyz.close()

	with open("atomtypes_tmp.dat","w") as atomtypes_tmp:
		for i in range(len(TmpAtoms)):
			atomtypes_tmp.write("{}-{}\n".format(AtomTypes[i],AtomicCharges[i]))
		atomtypes_tmp.close()

	with open("atomtypes_frozen.dat","w") as atomtypes_frozen:
		for i in range(len(TmpAtoms),len(FrozenAtoms)+len(TmpAtoms)):
			atomtypes_frozen.write("{}-{}\n".format(AtomTypes[i],AtomicCharges[i]))
		atomtypes_frozen.close()
