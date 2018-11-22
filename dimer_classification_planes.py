#!/usr/bin/env python3

import numpy as np
from sys import argv
import time
import cryspy.scripts.assign_charges as ac
import cryspy.io.read_file as rf
from sklearn.metrics.pairwise import pairwise_distances

start = time.time()

def open_check(xyzfile):
    """
    xyz file is checked to ensure that the number of coordinates matches the
    number of atoms.

    Parameters
    ----------
    xyzfile : The input xyz file to be read in

    Returns
    -------
    readfile : list of xyz file contents

    """

    readfile = open(xyzfile,"r").read().splitlines()

    number_atoms = int(readfile[0])
    length = len(readfile[2:])

    #Error message
    if number_atoms != length:
        raise ValueError("There are {0} atoms specified at the top of the file "
        "but {1} coordinates".format(number_atoms, length))

    return readfile

def distance(atom1, atom2):
    """Finds the distance between two atom coordinates"""
    vector = np.array(atom2)-np.array(atom1)
    #function finds euclidean distance between 2 atoms
    return np.linalg.norm(vector)

def write_xyz_kind(in_name, atoms):

    """
    Write an xyz file. (Works with function that assigns "kinds".)
    Needed for write_xyz_dimers to work

    Parameters
    ----------
    in_name : string
        Name of the xyz file. Include the file extension, e.g. "molecule.xyz"
    atoms : list of strings
        Atoms to write

    """
    out_file = open(in_name, "w")
    out_file.write(str(len(atoms)) + "\n")
    out_file.write(in_name + "\n")

    for atom in atoms:
        str_atom = str(atom)
        out_file.write(str_atom[:-8] +"\n")
    out_file.close()
    return

def write_xyz_dimers(unique_dimers):
    """Writes a separate xyz file for each dimer found."""
    for i,dimer in enumerate(unique_dimers):
        name=str(i)+'.xyz'
        flatten = [item for items in dimer for item in items]
        write_xyz_kind(name,flatten)

def less_than_threshold(A, B, threshold):
    for a,b in zip(A,B):
        if abs(a-b) > threshold: # distances are different
            return False
    return True

def find_connections(in_mono):
    """
    Assigns a 'kind' to each atom in an xyz file.

    Parameters
    ----------
    in_mono : atom coordinates as cryspy.utils.mol.Mol object.

    """
    first_connect = ac.detect_1_connect(in_mono, None, 1.8)
    total_connect = ac.complete_expand(first_connect)

    for i, atom in enumerate(in_mono):
        atom.set_connectivity(in_mono, total_connect[i])

def molecule_to_atoms(mol_object):
    """
    Takes a molecule object and returns the cartesian coordinates of each atom
    in the molecule.

    Parameters
    ----------
    mol_object : molecule object whose atom coordinates are to be found

    Returns
    -------
    molecule_positions : list of lists of np arrays of lists of coordinates
    Nmolecules x Natoms x array x coordinates lists
    """
    molecule_positions = []
    for molecule in mol_object:
            molecule_positions.append(molecule.get_pos())

    return molecule_positions

def dimers_centroid(molecules, cd):
    """
    Generates list of dimers based on centroid distances cd.
    (works with mol object)

    Parameters
    ----------
    cd : float
         Length between centroids of molecules

    Returns
    -------
    dimers : list of lists of strings
        2 x N list of [elem x y z] strings. They represent pairs of molecules

    """
    dimers=[] #list of dimers found
    for ref_number, ref_molecule in enumerate(molecules):
        for next_number, next_molecule in enumerate(molecules[ref_number+1:]):
            if not np.array_equal(molecule_to_atoms(ref_molecule), molecule_to_atoms(next_molecule)):
                centroid_ref=np.mean(molecule_to_atoms(ref_molecule),axis=0)
                centroid_next=np.mean(molecule_to_atoms(next_molecule),axis=0)
                if distance(centroid_ref, centroid_next)<= cd:
                    dimer=[ref_molecule,next_molecule]
                    dimers.append(dimer)
    return dimers

def inter_distances_dictionary(dimers):
    """
    Returns intermolecular distances for a dimer.

    Parameters
    ----------
    dimers : list of lists of mol object
    Ndimers x 2molecules x Natoms per molecule

    input_rounding_number : number of decimal places for distance result

    Returns
    -------
    dictionaries : list of dimers in a dictionary
    Ndimers x dictionary of atoms kind and interatomic distances

    """
    dictionaries=[]
    for dim_no, dim in enumerate(dimers):
        kinds_distances={}
        for atom_no_A, atom_A in enumerate(dim[0]):
            for atom_no_B, atom_B in enumerate(dim[1]):
                if atom_A.kind==atom_B.kind:
                    dist=distance(atom_A.get_pos(),atom_B.get_pos())
                    if atom_A.kind in kinds_distances.keys():
                        kinds_distances[atom_A.kind].append(dist)
                        kinds_distances[atom_A.kind]=sorted(kinds_distances[atom_A.kind])
                    else:
                        kinds_distances[atom_A.kind]=[dist]
        dictionaries.append(kinds_distances)
    return dictionaries

def purge(list_of_mols):
    """Remove incomplete molecules from list of Mols"""
    out_mols = []
    max_len = 0
    for mol in list_of_mols:
        if len(mol) > max_len:
            out_mols = []
            max_len = len(mol)
        if len(mol) == max_len:
            out_mols.append(mol)
    return out_mols

def SVD(X):
    """
    Defines a plane for a list of atom coordinates.
    Singular value decomposition method.
    Source: https://gist.github.com/lambdalisue/7201028
    """
    # Find the average of points (centroid) along the columns
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, Vt = np.linalg.svd(CX)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = Vt[-1]

    # Extract a, b, c, d coefficients.
    x0, y0, z0 = C
    a, b, c = N
    d = -(a * x0 + b * y0 + c * z0)

    return a, b, c, d

def project_to_plane(a,b,c,d,point):
    """
    Projects a point on to a plane of best fit.
    """
    m,n,o=point[0],point[1],point[2]
    t=(-d-a*m-b*n-c*o)/((a**2)+(b**2)+(c**2))
    projection=np.array([m+(t*a),n+(t*b),o+(t*c)])
    return projection

def squared_distance_to_plane(a,b,c,d,point):
    """
    Squared distance between a point and a plane.

    """
    x,y,z=point.T
    N=np.array([a,b,c])
    return np.mean(((a*x+b*y+c*z+d)**2/np.linalg.norm(N)))

def coordinate_matrix(atoms):
    coords=np.zeros((len(atoms),3))
    for i in range(len(atoms)):
        coords[i,0]=atoms[i].x
        coords[i,1]=atoms[i].y
        coords[i,2]=atoms[i].z
    return coords


def two_longest_axes(distance_matrix):
    """
    Finds the indices of the atoms that have the 2 longest pairwise distances in a distance matrix.

    argsort flattens and sorts the distance matrix.
    Unravel index then obtains the location of the index of the atoms that connect the longest distances
    in the matrix.
    The tuple is then returned as a transposed array such that the indices are paired
    if they are the longest or second longest distance within the distance matrix.

    Parameters
    ----------
    distance_matrix : NxN matrix of pairwise distances

    Returns
    -------
    an array of 2 arrays of 2 indices(int)

    """
    sort_index = np.argsort(distance_matrix, axis = None)
    return np.array(np.unravel_index(sort_index[[-2,-4]],distance_matrix.shape)).T

def midpoints(coordinates_1, coordinates_2):
    """
    Finds the midpoints of the shortest distances either side of a molecule.

    Parameters
    ----------
    coordinates_1 : 1x2 array of 1x3 coordinates The coordinates of the 2 points that contain the longest intramolecular distance.
    coordinates_2 : 1x2 array of 1x3 The coordinates of the 2 points that contain second longest intramolecular distance.

    Returns
    -------
    mid0 : The midpoint between the first coordinate of coordinates_1
    and the closest coordinate from coordinates_2.
    mid1 : The midpoint between the second coordinate of coordinates_1
    and the closest coordinate from coordinates_2.

    """
    c2min_0=coordinates_2[np.argmin([distance(coordinates_1[0],i) for i in coordinates_2])]
    c2min_1=coordinates_2[np.argmin([distance(coordinates_1[1],i) for i in coordinates_2])]
    mid0=np.mean([coordinates_1[0],c2min_0],axis=0)
    mid1=np.mean([coordinates_1[1],c2min_1],axis=0)
    return mid0, mid1

def point_to_line(line_1, line_2, point_project):
    """line_1 and line_2 are 2 points that define a line, point_project is a point to be projected onto the line"""
    A=point_project-line_1
    B=line_2-line_1
    projection = line_1 + np.dot(A,B)/np.dot(B,B)*B

    return projection

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """Angle between 2 vectors."""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

def slip_angle(centroid_1, centroid_2, centroid_perpendicular):
    """
    Finds the slip angle between 2 molecules.

    This is the angle between the centroid of the first molecule, the centroid
    of the second molecule and the vector perpendicular to the plane that contains
    the first molecule.
    """
    centroid_v = centroid_1-centroid_2
    slip = angle_between(centroid_perpendicular, centroid_v)

    if slip > 90:
        slip = 180-slip
    else:
        pass
    return slip


#finding unique dimers

read_file = argv[1]
#threshold = float(argv[2])
in_file = open_check(read_file)

clust = rf.mol_from_file(read_file)

segregated = clust.segregate()

list_of_molecules = purge(segregated)
for mol in list_of_molecules:
    find_connections(mol)
dimers = dimers_centroid(list_of_molecules, 10)

dictionaries=inter_distances_dictionary(dimers)

unique_dims = [dimers[0]]
unique_dictionaries = [[dictionaries[0], 1]]

# loop over all dimers
for i, dimer_dictionary in enumerate(dictionaries[1:]):
    unique = True
    for crosschecknumber,cross_check in enumerate(unique_dictionaries):
        for kind, dict_distance in dimer_dictionary.items():
            if less_than_threshold(dimer_dictionary[kind],cross_check[0][kind], threshold):
                unique = False
                cross_check[1]+=1
                break
        # if it's still unique after the checks
    if unique:
        unique_dims.append(dimers[i + 1])
        unique_dictionaries.append([dimer_dictionary, 1])

#Classifying Dimers

coordinate_matrices=[[coordinate_matrix(mon) for mon in dim] for dim in unique_dims]
distance_matrices=[[pairwise_distances(mon) for mon in dim] for dim in coordinate_matrices]

long_axes_angles=[]
short_axes_angles=[]
slip_angles=[]

for dim_index,dimer in enumerate(distance_matrices):
    long_axes_dimers=[]
    short_axes_dimers=[]
    centroids_dimer=[]
    perp_axes=[]

    for mon_index,monomer in enumerate(dimer):
        centroid_mon = np.mean(coordinate_matrices[dim_index][mon_index], axis=0)
        long1,long2 = two_longest_axes(monomer)
        coords1,coords2 = coordinate_matrices[dim_index][mon_index][long1],coordinate_matrices[dim_index][mon_index][long2]
        m1,m2=midpoints(coords1,coords2)
        a,b,c,d=SVD(coordinate_matrices[dim_index][mon_index])
        projected_m1=project_to_plane(a,b,c,d,m1)
        projected_m2=project_to_plane(a,b,c,d,m2)
        if distance(projected_m1,centroid_mon) < distance(projected_m2,centroid_mon):
            long_axis=projected_m2-projected_m1
        else:
            long_axis=projected_m1-projected_m2
        centroid_pro = point_to_line(projected_m1, projected_m2,centroid_mon)
        short_axis = centroid_pro - centroid_mon
        perp_axis = np.cross(short_axis, long_axis)


        perp_axes.append(perp_axis)
        centroids_dimer.append(centroid_mon)
        long_axes_dimers.append(long_axis)
        short_axes_dimers.append(short_axis)

    long_axes_angles.append(angle_between(long_axes_dimers[0],long_axes_dimers[1]))
    short_axes_angles.append(angle_between(short_axes_dimers[0],short_axes_dimers[1]))
    slip_angles.append(slip_angle(centroids_dimer[0],centroids_dimer[1], perp_axes[0]))

    perp_axes=[]
    centroids_dimer=[]
    long_axes_dimers = []
    short_axes_dimers = []


for index,(long_angle,short_angle,slip_angle) in enumerate(zip(long_axes_angles, short_axes_angles, slip_angles)):
    #print("Dimer {0:<5} Long:{1:<7.0f} Short:{2:<7.0f} Slip:{3:<7.0f}".format(index,long_angle,short_angle,slip_angle))
    print("{0:<7.0f} {1:<7.0f} {2:<7.0f}".format(long_angle, short_angle,slip_angle))



#write_xyz_dimers(unique_dims)



end = time.time()
#print("\nTotal time: {}s".format(round((end - start), 1)))
