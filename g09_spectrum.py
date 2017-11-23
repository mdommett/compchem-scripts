#! /usr/bin/env python

import numpy as np
import math
from sys import argv
from os import system
import argparse
#########
#########
"""
see http://gaussian.com/uvvisplot/ for implementation details

"""
parser = argparse.ArgumentParser()
parser.add_argument("input",help="Log file of Gaussian 09 TD job", type=str, nargs='*')
parser.add_argument("-gnu",help="Plot a spectrum using gnuplot",action="store_true")
parser.add_argument("-mpl",help="Plot a spectrum using matplotlib",action="store_true")
parser.add_argument("-sticks",help="Plot the stick spectrum",action="store_true")
parser.add_argument("-sd",help="Standard deviation (in eV)",default=0.4,type=float)
parser.add_argument("-r",help="Min and max values for the spectrum (in nm)",nargs=2,type=int)
parser.add_argument("-save",help="Save spectrum", type=str)
args=parser.parse_args()

def read_es(file):
    energies=[]
    os_strengths=[]
    for line in file:
        if " Excited State " in line:
            energies.append(float(line.split()[6]))
            os_strengths.append(float(line.split()[8][2:]))
    return energies,os_strengths

def abs_max(f,lam,ref):
    a=1.3062974e8
    b=f/(1e7/3099.6)
    c=np.exp(-(((1/ref-1/lam)/(1/(1240/args.sd)))**2))
    return a*b*c

def gnu_plot(xaxis,yaxis):
    with open("data","w") as d:
        for i in range(len(xaxis)):
            d.write("{0} {1} {2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

    with open("plot","w") as p:
        p.write("set xlabel \"Energy (nm)\"\nset ylabel \"Absorption Coeff.(E)\"\n")
        p.write("plot 'data' using 2:3 title '' with lines lt 1 lw 2,\\\n")
        p.close()

    system("gnuplot -persist plot")
    system("rm -f plot data")

    return

def mpl_plot(xaxis,yaxis):
    plt.scatter(xaxis,yaxis,s=2,c="r")
    plt.plot(xaxis,yaxis,color="k")
    plt.xlabel("Energy (nm)")
    plt.ylabel("$\epsilon$ (L mol$^{-1}$ cm$^{-1}$)")
    return

colours=["red","blue","green","orange"]
if args.mpl:
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

for n,f in enumerate(args.input):
    infile=open(f,"r").read().splitlines()
    energies,os_strengths=read_es(infile)

    if args.r:
        x=np.linspace(max(args.r),min(args.r),1000)

    else:
        x=np.linspace(max(energies)+200,min(energies)-200,1000)

    sum=[]
    for ref in x:
        tot=0
        for i in range(len(energies)):
            tot+=abs_max(os_strengths[i],energies[i],ref)
        sum.append(tot)

    stick_intensities=[abs_max(os_strengths[i],energies[i],energies[i]) for i in range(len(energies))]

    if args.gnu:
        gnu_plot(x,sum)

    elif args.mpl:
        plt.scatter(x,sum,s=2,c=colours[n])
        plt.plot(x,sum,color=colours[n],label=f[:-4])
        plt.xlabel("Energy (nm)")
        plt.ylabel("$\epsilon$ (L mol$^{-1}$ cm$^{-1}$)")
        if args.sticks:
            for i in range(len(energies)):
                plt.plot((energies[i],energies[i]),(0,stick_intensities[i]),colours[n])




    else:
        exit("Please specify -mpl or -gnu to plot with matplotlib or gnuplot")
if args.r:
    plt.xlim(min(args.r),max(args.r))
plt.legend()
if args.save:
    plt.savefig(args.save+".pdf")
plt.show()
