#!/usr/bin/env python
from cryspy.io import read_file as rf
from sys import argv
import numpy as np
from cryspy.io import edit_file as ef
from sklearn.metrics.pairwise import pairwise_distances

# usage
# dimer_classification.py dimer.xyz


def coordinate_matrix(atoms):
    coords=np.zeros((len(atoms),3))
    for i in range(len(atoms)):
        coords[i,0]=atoms[i].x
        coords[i,1]=atoms[i].y
        coords[i,2]=atoms[i].z
    return coords

def centroid(coordinates):
    length = coordinates.shape[0]
    sum_x = np.sum(coordinates[:, 0])
    sum_y = np.sum(coordinates[:, 1])
    sum_z = np.sum(coordinates[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])

def long_axis(distance_matrix):
    # Gets the maximum argument, which is a flattened index, and gives it back as a 2d index
    return np.array(np.unravel_index(np.argmax(distance_matrix, axis=None), distance_matrix.shape))

def distance(pointa,pointb):
    return np.sqrt(np.sum([(b-a)**2 for a,b in zip(point1,point2)]))


def align_z(coordinates,point1,point2):
    """first need to move a point to the origin
    then rotate into XZ plane
    """
    # translate chosen point to origin
    translated=coordinates-point1
    # long axis vector
    vector=point2-point1
    # unit_vector  of vector
    unit_vector=vector/np.sqrt(np.dot(vector, vector))
    # z-axis vector
    z=np.array([0,0,1])
    # cos(theta) of angle to between current vector and z-axis
    c=costheta(z,vector)
    rot_angle=(np.arccos(c))
    #print(np.degrees(rot_angle))
    #  rotation axis that connects vector and z-axis
    rot_axis=np.cross(z,vector)/magnitude(np.cross(z,vector))
    x,y,z=rot_axis
    C=1-c

    s=np.sin(rot_angle)
    # https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    # rotation matrix :
    rmat = np.array([[x*x*C+c,x*y*C-z*s,x*z*C+y*s],
                    [y*x*C+z*s,y*y*C+c,y*z*C-x*s],
                    [z*x*C-y*s,z*y*C+x*s,z*z*C+c]])
    new_coords=np.array([np.dot(i,rmat) for i in translated])
    return new_coords

def align_to_vector(coordinates,vector,ref_vector):
    """first need to move a point to the origin
    then rotate into XZ plane
    """

    unit_vector=vector/np.sqrt(np.dot(vector, vector))

    # cos(theta) of angle to between current vector and ref-vector
    c=costheta(ref_vector,vector)
    rot_angle=(np.arccos(c))

    #  rotation axis that connects vector and z-axis
    rot_axis=np.cross(ref_vector,vector)/magnitude(np.cross(ref_vector,vector))
    x,y,z=rot_axis
    C=1-c

    s=np.sin(rot_angle)
    # https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    # rotation matrix :
    rmat = np.array([[x*x*C+c,x*y*C-z*s,x*z*C+y*s],
                    [y*x*C+z*s,y*y*C+c,y*z*C-x*s],
                    [z*x*C-y*s,z*y*C+x*s,z*z*C+c]])
    new_coords=np.array([np.dot(i,rmat) for i in coordinates])
    return new_coords

def translate(coordinates,point):
    return coordinates-point

def magnitude(vector):
    return np.linalg.norm(vector)

def costheta(vector1,vector2):
    return np.dot(vector1,vector2)/(magnitude(vector1)*(magnitude(vector2)))

def distance(point1,point2):
    return np.sqrt(np.sum([(b - a)**2 for (a, b) in zip(point1,point2)]))

def dihedral_angle(p0,p1,p2,p3):
    # q0 = p1-p0 vector connecting p0 to p1
    q0 = p1 - p0
    q1 = p2 - p1
    q2 = p3 - p2

    q0q1=np.cross(q0,q1) # vector normal to plane p0p1p2
    q1q2=np.cross(q1,q2) # normal to plane p1p2p3

    n0=q0q1/magnitude(q0q1) # unit normal
    n1=q1q2/magnitude(q1q2)

    u0=n1 # orthogonal unit vectors
    u2=q1/magnitude(q1)
    u1=np.cross(u0,u2)

    cos_theta=np.dot(n0,u0)
    sin_theta=np.dot(n0,u1)
    theta=-np.arctan2(sin_theta,cos_theta)
    return np.degrees(theta)

def get_snap_vector((mon_0_0,mon_0_1),(mon_1_0,mon_1_1)):
    """
    Pass two coordinates per molecule (the extreme coordinates for the long axis)
    to generate the best vector to overlap the molecules
    the molecule with points 3 and 4 is translated
    """

    # intermolecular vectors between the extreme points of each molecule
    extreme_vectors=np.array([mon_1_0-mon_0_0,mon_1_0-mon_0_1,
    mon_1_1-mon_0_0,mon_1_1-mon_0_1])
    # translate each extreme point in mon_1 to each extreme point in molecule 0
    mon_1_translations=np.array([np.array((mon_1_0,mon_1_1)-vector) for vector in extreme_vectors])
    # new intermolecular extreme vectors
    new_extreme_vectors=np.array(mon_1_translations-(mon_0_0,mon_0_1))
    # best translation is the one that minimises the distances between the intermolecular extremes
    snap_vector=extreme_vectors[np.array([abs(magnitude(vector)) for vector in new_extreme_vectors]).argmin()]
    return snap_vector

def detect_element_indices(atoms,symbol):
    """
    get indexes for specified element from an array of atom objects
    """
    detected_indices=[index for index,atom in enumerate(atoms) if atom.elem==symbol]
    return detected_indices

def detect_nearest_element_distance(distance_matrix,ref_atom_index):
    """
    """
    return np.min(distance_matrix[ref_atom_index,:][np.nonzero(distance_matrix[ref_atom_index,:])])

def detect_CO_indexes(atoms,distance_matrix):
    """
    finds the indexes of the C and O of the hydroxyl group
    """
    # Find the oxgyen inddices
    oxygen_indices=detect_element_indices(atoms,'O')
    # Elements closet to oxygen
    oxygen_closest_distances=[detect_nearest_element_distance(distance_matrix,i) for i in oxygen_indices]
    # Find which oxygen is the C-O-H oxygen by checking if H is the nearest element to O
    for indice_index,distance in enumerate(oxygen_closest_distances):
        closest_atom_index=int(np.where(distance_matrix[oxygen_indices[indice_index],:]==distance)[0])
        closest_atom=atoms[closest_atom_index].elem
        if closest_atom=="H":
            #C=closest_atom_index
            O=oxygen_indices[indice_index]
    # Find C index by getting the second smallest distance value for COH (C-O) and placing it
    # at the third position (remember that zeros are included in distance)
    CO_distance=np.partition(distance_matrix[O,:],2)[2]
    # Get C index
    C=int(np.where(distance_matrix[O,:]==CO_distance)[0])
    return C,O

def detect_specific_nearest_element(distance_matrix,ref_atom_index,atom_indexes):
    distances=[distance_matrix[ref_atom_index,atom_index] for atom_index in atom_indexes]
    nearest_index=np.argmin(distances)
    other_index=np.argmax(distances)
    return atom_indexes[nearest_index],atom_indexes[other_index]


if __name__=='__main__':

    # load the dimer xyz coordinates
    dimer=rf.read_pos(argv[1])
    natoms_dimer=len(dimer)
    natoms_monomer=natoms_dimer/2
    # dimer must be arranged in sequence of monomers in xyz
    monomer_0=dimer[:natoms_monomer]
    monomer_1=dimer[natoms_monomer:]

    # get coordinates of each atom for each monomer
    mon_0=coordinate_matrix(dimer[0:natoms_monomer])
    mon_1=coordinate_matrix(dimer[natoms_monomer:])
    # get dimer coordinates
    dim=np.concatenate((mon_0,mon_1),axis=0)

    # distance matrices for each monomer
    mon_0_distances=pairwise_distances(mon_0)
    mon_1_distances=pairwise_distances(mon_1)

    # Get the atom indices of the C and O of the hydroxyl group
    # these are ordered by atom number
    mon_0_C,mon_0_O=(detect_CO_indexes(monomer_0,mon_0_distances))
    mon_1_C,mon_1_O=(detect_CO_indexes(monomer_1,mon_1_distances))
    mon_0_C_ab_distances=np.partition(mon_0_distances[mon_0_C],[2,3])[[2,3]]
    mon_1_C_ab_distances=np.partition(mon_1_distances[mon_1_C],[2,3])[[2,3]]
    mon_0_C_a,mon_0_C_b=[int(np.where(mon_0_distances[mon_0_C]==distance)[0]) for distance in mon_0_C_ab_distances]
    mon_1_C_a,mon_1_C_b=[int(np.where(mon_1_distances[mon_1_C]==distance)[0]) for distance in mon_1_C_ab_distances]
    mon_0_C_axis=mon_0[mon_0_C_b]-mon_0[mon_0_C_a]
    mon_1_C_axis=mon_1[mon_1_C_b]-mon_1[mon_1_C_a]
    # get the extreme atoms for each monomer, ordered by atom number
    mon_0_extreme_atoms_ordered=long_axis(mon_0_distances)
    mon_1_extreme_atoms_ordered=long_axis(mon_1_distances)#+natoms_monomer # back to index of original dimer

    # reorder these so that the first extreme atom is always closest to the carbonyl and the second is always
    # furthest from carbonyl
    mon_0_extreme_atoms=detect_specific_nearest_element(mon_0_distances,mon_0_C,mon_0_extreme_atoms_ordered)
    mon_1_extreme_atoms=detect_specific_nearest_element(mon_1_distances,mon_1_C,mon_1_extreme_atoms_ordered)

    # get the coordinates of each extreme atom index
    mon_0_extreme_coords=np.array([mon_0[mon_0_extreme_atoms[0]],mon_0[mon_0_extreme_atoms[1]]])
    mon_1_extreme_coords=np.array([mon_1[mon_1_extreme_atoms[0]],mon_1[mon_1_extreme_atoms[1]]])

    # define the long axis vector for each monomer
    mon_0_long_axis=mon_0_extreme_coords[1]-mon_0_extreme_coords[0]
    mon_1_long_axis=mon_1_extreme_coords[1]-mon_1_extreme_coords[0]

    # define the short axis vector for each monomer
    mon_0_CO_axis=mon_0[mon_0_O]-mon_0[mon_0_C]
    mon_1_CO_axis=mon_1[mon_1_O]-mon_1[mon_1_C]

    # calculate the long axis angle
    long_axis_angle=np.degrees(np.arccos(costheta(mon_0_long_axis,mon_1_long_axis)))

    # calculate the short axis angle
    CO_angle=np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_CO_axis)))
    # calculate the angle between long and short (not used in analysis)
    CO_long_axis_angle=np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_long_axis)))
    #print(np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_long_axis))))

    # align to CO to Z-axis
    aligned_dim=align_z(dim,mon_0[mon_0_C],mon_0[mon_0_O])
    aligned_mon_0=aligned_dim[:natoms_monomer]
    aligned_mon_1=aligned_dim[natoms_monomer:]
    mon_0_C_axis=aligned_mon_0[mon_0_C_b]-aligned_mon_0[mon_0_C_a]
    mon_0_CO_axis=aligned_mon_0[mon_0_O]-aligned_mon_0[mon_0_C]
    mon_0_perpendicular_axis=np.cross(mon_0_C_axis,mon_0_CO_axis)
    # align to C axis to Y-axis
    aligned_dim=align_to_vector(aligned_dim,mon_0_C_axis,np.array([0,1,0]))
    aligned_mon_0=aligned_dim[:natoms_monomer]
    aligned_mon_1=aligned_dim[natoms_monomer:]

    mon_0_C_axis=aligned_mon_0[mon_0_C_b]-aligned_mon_0[mon_0_C_a]
    mon_1_C_axis=aligned_mon_1[mon_1_C_b]-aligned_mon_1[mon_1_C_a]
    C_axis_angle=np.degrees(np.arccos(costheta(mon_0_C_axis,mon_1_C_axis)))

    mon_0_CO_axis=aligned_mon_0[mon_0_O]-aligned_mon_0[mon_0_C]
    mon_1_CO_axis=aligned_mon_1[mon_1_O]-aligned_mon_1[mon_1_C]
    CO_angle=np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_CO_axis)))
    perp_vector=np.cross(mon_0_C_axis,mon_0_CO_axis)

    slip_angle=np.degrees(np.arccos(costheta(perp_vector,aligned_mon_1[mon_1_C]-aligned_mon_0[mon_0_C])))
    mon_0_centroid=centroid(aligned_mon_0)
    mon_1_centroid=centroid(aligned_mon_1)
    print(mon_0_centroid)
    print(mon_1_centroid)
    y_slip=abs(mon_1_centroid[1]-mon_0_centroid[1])
    aligned_mon_0_long_axis=aligned_mon_0[mon_0_extreme_atoms[1]]-aligned_mon_0[mon_0_extreme_atoms[0]]
    mon_0_y_length=abs(aligned_mon_0_long_axis[1])
    mon_0_centroid_axis=(mon_0_centroid+np.array([1,0,0]))-mon_0_centroid
    centroid_angle=np.degrees(np.arccos(costheta(mon_0_centroid_axis,mon_1_centroid-mon_0_centroid)))
    slip_angle_centroid=np.degrees(np.arccos(costheta(perp_vector,mon_1_centroid-aligned_mon_0[mon_0_C])))
    if slip_angle >90:
        slip_angle=180-slip_angle
    if slip_angle_centroid >90:
        slip_angle_centroid=180-slip_angle_centroid
    if centroid_angle >90:
        centroid_angle=180-centroid_angle
    print("{0:>7.3f} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4:>7.3f} {5:>7.3f} {6:>7.3f}".format(C_axis_angle,CO_angle,slip_angle,slip_angle_centroid,centroid_angle,y_slip,mon_0_y_length))
    for i,j in enumerate(dimer):
        j.x,j.y,j.z=aligned_dim[i]
    ef.write_xyz("{}_aligned.xyz".format(argv[1][:-4]),dimer)
    ############################################################################################
    """
    ## DISTANCE ANALYSIS NOT USED/NEEDED CURRENTLY##
    # snap vector aligns the long axes of each monomer
    snap_vector=get_snap_vector(mon_0_extreme_coords,mon_1_extreme_coords)
    magnitude_snap_vector=magnitude(snap_vector)
    # align vectors
    mon_1_snapped=mon_1-snap_vector
    snapped_dim=np.concatenate((mon_0,mon_1_snapped))
    translated_dim=translate(dim,mon_0_extreme_coords[0])


    # Align long axis of mon_0 along z-axis, so that z-axis slip  and centroid z-slip
    #can easily be identified

    aligned_dim=align_z(dim,mon_0_extreme_coords[0],mon_0_extreme_coords[1])
    aligned_mon_0=aligned_dim[:natoms_monomer]
    aligned_mon_1=aligned_dim[natoms_monomer:]


    mon_0_centroid=centroid(aligned_mon_0)
    mon_1_centroid=centroid(aligned_mon_1)
    # z-slip of centroid
    centroid_z_slip=mon_1_centroid[2]-mon_0_centroid[2]
    # centroid distance
    centroid_distance=np.linalg.norm((mon_1_centroid-mon_0_centroid))
    mon_z_length=abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    centroid_z_slip_norm=centroid_z_slip/mon_z_length
    CO_slip=aligned_mon_1[mon_1_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[0]][2]
    CO_N_slip=aligned_mon_1[mon_1_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2]
    N_slip=aligned_mon_1[mon_1_extreme_atoms[1]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2]
    CO_slip_norm=CO_slip/mon_z_length
    CO_N_slip_norm=CO_N_slip/abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    N_slip_norm=N_slip/abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    #print("{0:>7.3f} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4:>7.3f} {5:>7.3f} {6:>7.3f} {7:>7.3f} {8:>7.3f}".format(CO_slip,CO_N_slip,N_slip,centroid_distance,centroid_z_slip,long_axis_angle,CO_angle,CO_long_axis_angle,mon_z_length))
    """
