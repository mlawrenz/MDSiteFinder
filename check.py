import numpy
import scipy.spatial as sp


def make():
    array1=numpy.array([1,1,1])
    array2=numpy.array([2,2,2])
    ref=numpy.array([4,4,4])
    new_matrix=numpy.zeros((2,2,3))
    new_matrix[0][0]=array1
    new_matrix[0][1]=array2
    new_matrix[1][0]=array1
    new_matrix[1][1]=array2
    ref_matrix=numpy.zeros((2,2,3))
    ref_matrix[0][0]=ref
    ref_matrix[0][1]=ref
    ref_matrix[1][0]=ref
    ref_matrix[1][1]=ref
    return new_matrix, ref_matrix

#output=sp.distance.euclidean(ref_matrix, new_matrix)
