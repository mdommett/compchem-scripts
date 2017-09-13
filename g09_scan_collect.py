#!/usr/bin/env python

from sys import argv,exit
import argparse
from periodic import element
"""
Parse the output of a Gaussian 09 relaxed geometry scan job

Will automatically save the corresponding energies, with options to
save the optimised geometries [-xyz] and plot a graph of the s0 and s1 states
[-p]

Michael Dommett
August 2017
m.dommett@qmul.ac.uk
"""


parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True,help="Log file of Gaussian 09 scan")
parser.add_argument("-xyz","--xyz",help="Save the optimised geometries in xyz format",action="store_true")
parser.add_argument('-p',"--plot",help="Plot a graph of the energies",action="store_true")
args=parser.parse_args()

def read_energies(infile):

    """
    Parses the ground and excited state energy for each optimisation step

    Accepts a G09 output file as an input

    Parameters
    ----------
    infile : list of str
        G09 logfile
    Returns
    -------
    scf_final : list of ground state energy values, in eV, relative to the
    initial ground state value
    s1_final : list of excited state energy values, in eV, relative to the
    initial ground state value
    """
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
    """
    Parses the logfile for the optimised geometry at each step, writes an output
    file with the geometries in xyz format

    Accepts a G09 output file as an input

    Parameters
    ----------
    infile : list of str
        G09 logfile
    """
    print "Saving each optimised geometry to {} ...\n".format(args.input+"-opt.xyz")
    fgeom=open(args.input+"-opt.xyz","w")
    for i,j in enumerate(infile):
        if j.split()[0]=="NAtoms=":
            natoms=int(j.split()[3])
            break
    for i,j in enumerate(infile):
        if j == " Number     Number       Type             X           Y           Z":
            hold_geom=[]
            count = int(i)
            hold_geom.append(lines[count+2:count+natoms+2])
        if j==" Optimization completed.":
            fgeom.write("{}\n\n".format(natoms))
            for geom in hold_geom:
                for line in geom:
                    splt = line.split()
                    symbol,x,y,z=(element(splt[1]).symbol,float(splt[3]),float(splt[4]),float(splt[5]))
                    fgeom.write("{0:<2} {1:>13.9f} {2:>13.9f} {3:>13.9f}\n".format(symbol,x,y,z))
    return

def plot_energies(scf_final,s1_final):
    """
    Plots the ground and excited states energies which are generated by
    the read_energies function

    Parameters
    ----------
    scf_final : list of ground state energy values, in eV, relative to the
    initial ground state value
    s1_final : list of excited state energy values, in eV, relative to the
    initial ground state value
    Returns
    -------
    Plot of ground and excited state energy vs scan step number
    """
    fig,ax = plt.subplots()
    x=range(len(s1_final))
    ax.plot(x,scf_final,'o-b',linewidth=2)
    ax.plot(x,s1_final,'o-r',linewidth=2)
    plt.xlabel('Step Number', fontsize=12)
    plt.ylabel('Energy (eV)', fontsize=12)
    #plt.show()
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
    plot_energies(scf_final,s1_final)
    plt.show()
