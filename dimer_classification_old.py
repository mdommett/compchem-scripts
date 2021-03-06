#!/usr/bin/env python
import read_file as rf
from sys import argv
import numpy as np
import edit_file as ef
from sklearn.metrics.pairwise import pairwise_distances
#np.seterr(divide='ignore', invalid='ignore') # divide by zero


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

def rotation_matrix(axis, theta):

    """Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians."""

    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))

    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotate_x(v,theta):
    rot=np.array([[1,0,0],
                [0,np.cos(theta),-np.sin(theta)],
                [0,np.sin(theta),np.cos(theta)]])
    return np.dot(rot,v)

def rotate(axis,theta):
    axis = axis/np.sqrt(np.dot(axis, axis))
    a11=np.cos(theta)+(axis[0]**2*(1-np.cos(theta)))
    a12=axis[0]*axis[1]*(1-np.cos(theta))-axis[2]*np.sin(theta)
    a13=axis[0]*axis[2]*(1-np.cos(theta))+axis[1]*np.sin(theta)
    a21=axis[1]*axis[0]*(1-np.cos(theta))+axis[2]*np.sin(theta)
    a22=np.cos(theta)+(axis[1]**2*(1-np.cos(theta)))
    a23=(axis[1]*axis[2]*(1-np.cos(theta)))-axis[0]*np.sin(theta)
    a31=(axis[2]*axis[0]*(1-np.cos(theta)))-axis[1]*np.sin(theta)
    a32=(axis[2]*axis[1]*(1-np.cos(theta)))+axis[0]*np.sin(theta)
    a33=np.cos(theta)+(axis[2]**2*(1-np.cos(theta)))
    R=np.array([[a11,a12,a13],
                [a21,a22,a23],
                [a31,a32,a33]])
    return R

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

    dimer=rf.read_pos(argv[1])
    natoms_dimer=len(dimer)
    natoms_monomer=natoms_dimer/2
    monomer_0=dimer[:natoms_monomer]
    monomer_1=dimer[natoms_monomer:]

    mon_0=coordinate_matrix(dimer[0:natoms_monomer])
    mon_1=coordinate_matrix(dimer[natoms_monomer:])
    dim=np.concatenate((mon_0,mon_1),axis=0)

    mon_0_distances=pairwise_distances(mon_0)
    mon_1_distances=pairwise_distances(mon_1)
    mon_0_C,mon_0_O=(detect_CO_indexes(monomer_0,mon_0_distances))
    mon_1_C,mon_1_O=(detect_CO_indexes(monomer_1,mon_1_distances))

    mon_0_extreme_atoms_ordered=long_axis(mon_0_distances)
    mon_0_extreme_atoms=detect_specific_nearest_element(mon_0_distances,mon_0_C,mon_0_extreme_atoms_ordered)
    mon_1_extreme_atoms_ordered=long_axis(mon_1_distances)#+natoms_monomer # back to index of original dimer
    mon_1_extreme_atoms=detect_specific_nearest_element(mon_1_distances,mon_1_C,mon_1_extreme_atoms_ordered)
    mon_0_extreme_coords=np.array([mon_0[mon_0_extreme_atoms[0]],mon_0[mon_0_extreme_atoms[1]]])
    mon_1_extreme_coords=np.array([mon_1[mon_1_extreme_atoms[0]],mon_1[mon_1_extreme_atoms[1]]])
    print(mon_0_extreme_atoms,mon_1_extreme_atoms)
    mon_0_long_axis=mon_0_extreme_coords[1]-mon_0_extreme_coords[0]
    mon_0_CO_axis=mon_0[mon_0_O]-mon_0[mon_0_C]
    mon_1_long_axis=mon_1_extreme_coords[1]-mon_1_extreme_coords[0]
    mon_1_CO_axis=mon_1[mon_1_O]-mon_1[mon_1_C]
    #print(dihedral_angle(mon_0[mon_0_C],mon_0[mon_0_O],mon_1[mon_1_C],mon_1[mon_1_O]))
    #print(np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_0_long_axis))))
    long_axis_angle=np.degrees(np.arccos(costheta(mon_0_long_axis,mon_1_long_axis)))
    CO_angle=np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_CO_axis)))
    CO_long_axis_angle=np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_long_axis)))
    #print(np.degrees(np.arccos(costheta(mon_0_CO_axis,mon_1_long_axis))))
    #print(np.degrees(np.arccos(costheta(mon_0_long_axis,mon_1_long_axis))))

    #print(dihedral_angle(mon_0[mon_0_extremes[0]],mon_0[mon_0_extremes[1]],mon_1[mon_1_extremes[0]],mon_1[mon_1_extremes[1]]))


    snap_vector=get_snap_vector(mon_0_extreme_coords,mon_1_extreme_coords)
    magnitude_snap_vector=magnitude(snap_vector)
    mon_1_snapped=mon_1-snap_vector
    snapped_dim=np.concatenate((mon_0,mon_1_snapped))
    translated_dim=translate(dim,mon_0_extreme_coords[0])
    aligned_dim=align_z(dim,mon_0_extreme_coords[0],mon_0_extreme_coords[1])
    aligned_mon_0=aligned_dim[:natoms_monomer]
    aligned_mon_1=aligned_dim[natoms_monomer:]
    mon_0_centroid=centroid(mon_0)
    mon_1_centroid=centroid(mon_1)
    centroid_z_slip=mon_1_centroid[2]-mon_0_centroid[2]
    centroid_distance=np.linalg.norm((mon_1_centroid-mon_0_centroid))
    mon_z_length=abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    centroid_z_slip_norm=centroid_z_slip/mon_z_length
    CO_slip=aligned_mon_1[mon_1_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[0]][2]
    CO_N_slip=aligned_mon_1[mon_1_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2]
    N_slip=aligned_mon_1[mon_1_extreme_atoms[1]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2]
    CO_slip_norm=CO_slip/mon_z_length
    CO_N_slip_norm=CO_N_slip/abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    N_slip_norm=N_slip/abs(aligned_mon_0[mon_0_extreme_atoms[0]][2]-aligned_mon_0[mon_0_extreme_atoms[1]][2])
    print("{0:>7.3f} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4:>7.3f} {5:>7.3f} {6:>7.3f} {7:>7.3f} {8:>7.3f}".format(CO_slip,CO_N_slip,N_slip,centroid_distance,centroid_z_slip,long_axis_angle,CO_angle,CO_long_axis_angle,mon_z_length))

    #mon_0_new_long_axis=snapped_and_aligned_dim[mon_0_extreme_atoms[1]]-snapped_and_aligned_dim[mon_0_extreme_atoms[0]]
    #mon_1_new_long_axis=snapped_and_aligned_dim[mon_1_extreme_atoms[1]]-snapped_and_aligned_dim[mon_1_extreme_atoms[0]]
    #COangle=np.degrees(np.arccos(costheta(CO1,CO2)))
"""
    for i ,j in enumerate(dimer):
        j.x,j.y,j.z=snapped_and_aligned_dim[i]
    ef.write_xyz("{}_snapped_and_aligned.xyz".format(argv[1][:-4]),dimer)
    for i ,j in enumerate(dimer):
        j.x,j.y,j.z=aligned_dim[i]
    ef.write_xyz("{}_aligned.xyz".format(argv[1][:-4]),dimer)

    #print(z_CO)

    #aligned_dim=align_z(dim,mon_0_extreme_coords[0],mon_0_extreme_coords[1])

    #mon_0_new_long_axis=translated_dim[mon_0_extreme_atoms[1]]-translated_dim[mon_0_extreme_atoms[0]]
    #print(np.degrees(costheta(mon_0_new_long_axis,(1,0,0))))
    ###
    # dimer angle (pyr)
    cvec=translated_dim[37]-mon_0_new_long_axis[0]
    dvec=translated_dim[37]-mon_0_new_long_axis[1]
    pyr_angle=np.degrees(np.arccos(np.dot(np.cross(cvec,dvec),CO1)/np.dot(magnitude(np.cross(cvec,dvec)),magnitude(CO1))))
    ###


    ###
    #print(mon_1_snapped[mon_1_extreme_atoms[0]])
    #snapped_dim=np.concatenate((mon_0,mon_1_snapped))
    for i ,j in enumerate(dimer):
        j.x,j.y,j.z=aligned_dim[i]
    ef.write_xyz("{}_aligned.xyz".format(argv[1][:-4]),dimer)

    for no,vector in enumerate(extreme_vectors):
        #print(magnitude(vector))
        #theta = np.radians(angle)
        #rotated=np.array([np.dot(rotate(np.array([0,0,1]),theta), i) for i in new_coords])
        mon_1_snapped=mon_1-vector
        print(np.mean(pairwise_distances(mon_0,mon_1_snapped)))
        snapped_dim=np.concatenate((mon_0,mon_1_snapped))
        #print(distance(snapped_dim[19],snapped_dim[7]))
        #translated_dim=align_z(snapped_dim,mon_0_extreme_coords[0],mon_0_extreme_coords[1])
        for i,j in enumerate(dimer):
            j.x,j.y,j.z=snapped_dim[i]
        ef.write_xyz("{}.xyz".format(no),dimer)
        """
