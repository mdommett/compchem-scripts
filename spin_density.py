#!/usr/bin/env python

# Script to read and then  plot the spin density from a g09 output.

import numpy as np
from sys import argv,exit
from periodic import element
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",help="[Mull] Mulliken charges or [NBO] NBO charges", default="Mull")
parser.add_argument("-i","--input",required=True,help="Input file")
parser.add_argument("-p","--plot",help="Plot output",action='store_true')
args=parser.parse_args()

if args.type=="NBO":
    NBO_in = open(args.input,"r").read().splitlines()

    start = []
    stop = []

    count=0
    for i,j in enumerate(NBO_in):

        if j == "    Atom  No    Charge         Core      Valence    Rydberg      Total":

            start.append(i)
        if j == " =======================================================================":
            stop.append(i)
    alphasymbol,alphanumber,alphanatural_charge = [],[],[]
    betasymbol,betanumber,betanatural_charge = [],[],[]

    if len(start)==len(stop):
        for line in NBO_in[start[1]+2:stop[1]]:
            cols = line.split()
            alphasymbol.append(cols[0])
            alphanumber.append(int(cols[1]))
            alphanatural_charge.append(float(cols[2]))
        for line in NBO_in[start[2]+2:stop[2]]:
            cols = line.split()
            betasymbol.append(cols[0])
            betanumber.append(int(cols[1]))
            betanatural_charge.append(float(cols[2]))
    else:
        sys.exit("Error reading charges")
    density=[]
    if alphasymbol==betasymbol and alphanumber==betanumber:
        for i in range(len(alphasymbol)):
            print "{0:<2}  {1:>3}  {2:>8.5f}".format(alphasymbol[i], alphanumber[i], betanatural_charge[i]-alphanatural_charge[i])
            density.append(betanatural_charge[i]-alphanatural_charge[i])
    else:
        sys.exit("Error")

    if args.plot:
        fig,ax=plt.subplots()
        ax.set_title(args.input[:-4])
        ax.set_xlabel("Atom Number")
        ax.set_ylabel("Spin density")
        ax.bar(alphanumber,density)
        plt.savefig(str(args.input)[:-4])+"_spindensity.pdf"
        plt.show()
else:
    infile = open(args.input,"r").read().splitlines()
    for linenumber,line in enumerate(infile):
        if line==" Mulliken charges and spin densities:":
            start=linenumber+2
        if " Sum of Mulliken charges =" in line:
            stop=linenumber
            break
    spindensities=[float(i.split()[3]) for i in infile[start:stop]]

    if args.plot:
        fig,ax=plt.subplots()
        ax.set_title(args.input[:-4])
        ax.set_xlabel("Atom Number")
        ax.set_ylabel("Spin density")
        ax.bar(range(len(spindensities)),spindensities)
        plt.savefig(str(args.input[:-4])+"_spindensity.pdf")
        plt.show()
