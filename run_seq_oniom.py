#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv
import sys

finp=open("gaussian.inp","r").read().splitlines()
fcom=open("gaussian.com","w")
fout=open("gaussian.log","w")
tmpout=open("tmp.out","w")
tmpcom=open("tmp.com","r").read().splitlines()

print "Starting program\n"
tmpcom=open("tmp.com","r").read().splitlines()
nat=int(tmpcom[0])
print "Number of atoms: {}".format(nat)

with open("method.inp","r") as mthd:
    prin

