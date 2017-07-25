#!/usr/bin/env python

import numpy as np
from itertools import chain
import sys
from sys import argv

input_file = open(argv[1] + ".out","r")
infile = input_file.read().splitlines()
outfile = open(argv[1] + ".xyz","w")

for i,j in enumerate(infile):
    if j ==" Symbolic Z-Matrix:":
        count = int(i)
        
for i,j in enumerate(infile[count+1:]):
    if j =="":
        natoms =  int(i)-1
        break
        
for i,j in enumerate(infile):

    if j == "                           !   Optimized Parameters   !":
        count = int(i)
        break
        
for i,j in enumerate(infile[count:]): 
    if j == " Number     Number       Type             X           Y           Z":
        test = int(i)
        break