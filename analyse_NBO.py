#!/usr/bin/env python

from sys import argv,exit
import numpy as np

infiles0 = open(argv[1],"r")
infiles1 = open(argv[2],"r")
NBO_s0 = infiles0.read().splitlines()
NBO_s1 = infiles1.read().splitlines()

files = NBO_s0,NBO_s1
search_string = " Atom  No    Charge         Core      Valence    Rydberg      Total"

def get_line_number(list,int):
    
    for line in list:
        if search_string in line:
            break
        else:
            int+=1
    return 

current_line=1
get_line_number(NBO_s0,current_line)
print current_line


              
