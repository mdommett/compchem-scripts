#!/usr/bin/env python

"""
Plot one or more cross-section.dat NX file

Supports saving the file, overlapping any number of graphs in nm or eV
and Gaussian 09 oscillator strengths

31-07-2017 Miguel Rivera
"""

import sys
import matplotlib.pyplot as plt
from matplotlib import collections as matcoll
from cycler import cycler
import argparse
import numpy

def read_NX(in_name, in_unit):
    """
    Read an NX spectrum

    Pass it a list of strings ['file in_name','in_unit']

    Parameters
    ----------
    in_name : str
        The name of the file to read
    in_unit : str
        Unit of measurement, "eV" or "nm"
    Returns
    -------
    x : list of floats
        Energies (wavelengths) of the spectrum
    y : list of floats
        Intensities in arbitrary units
    """
    unit_n = 0
    if in_unit.lower() == "ev":
        unit_n = 0
    if in_unit.lower() == "nm":
        unit_n = 1
    with open(in_name) as data:
        lines = data.readlines()

    y = []
    x = []
    for line in lines[1:]:
        x.append(float(line.split()[unit_n]))
        y.append(float(line.split()[2]))
    y_max=max(y)
    y = [i/y_max for i in y]

    return x, y


def read_ex_gaussian(in_name):
    """
    Read a Gaussian .log file with excitations

    Returns the excitation energies and oscillator strengths

    Parameters
    ----------
    in_name : str
        Name of the Gaussian log file including extension
    Returns
    -------
    energies_ev = list of floats
        Excitation energies in ev
    energies_nm = list of floats
        Excitation energies in ev
    stren = list of floats
        Oscillator strengths

    """
    with open(in_name) as data:
        lines = data.readlines()

    energies_ev = []
    energies_nm = []

    stren = []

    for line in lines:
        if "Excited State" in line:
            energies_ev.append(float(line.split()[4]))
            energies_nm.append(float(line.split()[6]))
            osci = float(line.split()[8][2:])
            stren.append(osci)
    return energies_ev, energies_nm, stren


def read_exp(in_name):
    """
    Read an experimental data file in 2 column format

    Only in nm

    Parameters
    ----------
    in_name : str
        Name of the 2 column file to read
    Returns
    -------
    x : list of floats
        Energies in nm
    y : list of floats
        Intensities in arbitrary units

    """
    with open(in_name) as data:
        lines = data.readlines()

    x=[]
    y=[]
    for line in lines:
        if line.strip():
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))
    y_max=max(y)
    y = [i/y_max for i in y]

    return x, y

def main(unit,gaussian,experimental,pdf,jpg,show,nx_filename):
    # sets of spectrum data
    x_list = []
    y_list = []

    for file_name in nx_filename:
        x_temp, y_temp = read_NX(file_name, unit)
        x_list.append(x_temp)
        y_list.append(y_temp)
    fig, total_plot = plt.subplots()

    total_plot.set_prop_cycle(cycler('color', ['c', 'm', 'y', 'k']))
    for x_item, y_item in zip(x_list, y_list):
        total_plot.plot(x_item, y_item)

    # labels
    plt.ylabel('Cross section', fontsize=12)
    plt.xlabel('Energy in ' + unit, fontsize=12)
    # some options
    if gaussian:
        for g_set in gaussian:

            # convoluted bit to draw vertical lines
            lines = []
            e_ev, e_nm, osc = read_ex_gaussian(g_set)
            if unit.lower() == "ev":
                exci = e_ev
            if unit.lower() == "nm":
                exci = e_nm

            for i in range(len(exci)):
                pair = [(exci[i], 0), (exci[i], osc[i])]
                lines.append(pair)
            linecoll = matcoll.LineCollection(lines,colors=numpy.random.rand(3,))
            # for the second y axis
            g_plot = total_plot.twinx()
            g_plot.add_collection(linecoll)
    if experimental:
        for e_set in experimental:
            e_x, e_y = read_exp(e_set)
#            e_plot = total_plot.twinx()
            total_plot.plot(e_x, e_y)
    if pdf:
        plt.savefig("out.pdf")
    if jpg:
        plt.savefig("out.jpg")
    if show:
        plt.show()


if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u", "--unit", help="Pick the unit, eV or nm", default="nm", type=str)
    parser.add_argument("-g", "--gaussian", help="Gaussian excitation .log files",
                        default=[], type=str, nargs='*')
    parser.add_argument("-e", "--experimental", help="Table of two columns: nm and intensity",
                        default=[], type=str, nargs='*')
    parser.add_argument("-p", "--pdf", help="Save a out.pdf output",
                        action="store_true", default=False)
    parser.add_argument("-j", "--jpg", help="Save a out.jpg output",
                        action="store_true", default=False)
    parser.add_argument(
        "-s", "--show", help="Show the figure or not", action="store_true", default=True)
    parser.add_argument("nx_filename", help="Newton-X files to open",
                        default=["cross-section.dat"], type=str, nargs='*')

    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b','y'])))
    read_unit = args.unit
    if read_unit.lower() != "ev" and read_unit.lower() != "nm":
        raise ValueError("Invalid unit, please choose eV or nm")
    read_gaussian = args.gaussian
    read_experimental = args.experimental
    read_pdf = args.pdf
    read_jpg = args.jpg
    read_show = args.show
    read_nx_filename = args.nx_filename

    #call main
    main(read_unit,read_gaussian,read_experimental,read_pdf,read_jpg,read_show,read_nx_filename)
