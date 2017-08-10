#!/usr/bin/env python

from sys import argv,exit
import argparse
from periodic import element

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True,help="Log file of Gaussian 09 scan")
parser.add_argument("-xyz","--xyz",help="Save the optimised geometries in xyz format",action="store_true")
parser.add_argument('-p',"--plot",help="Plot a graph of the energies",action="store_true")
args=parser.parse_args()

def read_energies(infile):

    hold_scf=[]
    hold_s1=[]
    for line in infile:

        if  " SCF Done:" in line:
            scf=float(line.split()[4])

        if " Total Energy" in line:
            s1=float(line.split()[4])

        if line==" Optimization completed.":
            hold_scf.append(scf)
            hold_s1.append(s1)
            del scf,s1

    scf_final=[(i-hold_scf[0])*27.2114 for i in hold_scf]
    s1_final=[(i-hold_scf[0])*27.2114 for i in hold_s1]

    return scf_final,s1_final

def save_geometries(infile):
    print "Saving each optimised geometry to {} ...\n".format(args.input+"-opt.xyz")
    fgeom=open(args.input+"-opt.xyz","w")
    hold_geom=[]
    for i,j in enumerate(infile):
        if j.split()[0]=="NAtoms=":
            natoms=int(j.split()[3])
            break
    for i,j in enumerate(infile):
        if j == " Number     Number       Type             X           Y           Z":
            hold_geom=[]
            count = int(i)
            for line in lines[count+2:count+natoms+2]:
                splt = line.split()
                coords=(element(splt[1]).symbol,float(splt[3]),float(splt[4]),float(splt[5]))
                hold_geom.append(coords)
        if j==" Optimization completed.":
            fgeom.write("{}\n\n".format(natoms))
            for line in hold_geom:
                symbol,x,y,z=line
                fgeom.write("{0:<2} {1:>13.9f} {2:>13.9f} {3:>13.9f}\n".format(symbol,x,y,z))
    return

if __name__=='__main__':

    finp = open(args.input,"r").read().splitlines()
    fen = open(args.input+"-energies","w")
    lines = filter(None, (line.rstrip() for line in finp))

    scf_final,s1_final=read_energies(lines)
    print "Saving energies to {} ...\n".format(args.input+"-energies")
    if len(scf_final)==len(s1_final):
        for i in range(len(scf_final)):
            fen.write("{0:>2.3f}  {1:>2.3f}\n".format(scf_final[i],s1_final[i]))
    else:
        sys.exit("\n\nCheck the G09 output file, there is a problem\n\n")

    if args.xyz:
        save_geometries(lines)

if args.plot:
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots()
    x=range(len(s1_final))
    ax.plot(x,scf_final,'o-b',linewidth=2)
    ax.plot(x,s1_final,'o-r',linewidth=2)
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Energy (eV)', fontsize=12)
    plt.show()
