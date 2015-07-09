import sys
import glob
import os
import sys
import optparse
import pickle
import numpy

def assign(matrix, database, cutoff):
    distances=-1*numpy.ones(len(database))
    assignments=-1*numpy.ones(len(database), dtype=int)
    for j in xrange(matrix.shape[0]):
        d=matrix[j,:]
        ind=numpy.argmin(d)
        if not d[ind] < cutoff:
            pass
        else:
            assignments[j] = int(numpy.argmin(d))
            distances[j] = d[assignments[j]]
    return assignments, distances
 
def cluster(matrix, distance_cutoff, cluster_cutoff=None):
    if cluster_cutoff is None and distance_cutoff is None:
        raise ValueError("I need some cutoff criterion! both k and distance_cutoff can't both be none")
    if cluster_cutoff is None and distance_cutoff <= 0:
        raise ValueError("With k=None you need to supply a legit distance_cutoff")
    if distance_cutoff is None:
        # set it below anything that can ever be reached
        distance_cutoff = -1
    if cluster_cutoff is None:
        # set k to be the highest 32bit integer
        cluster_cutoff = sys.maxint

    distance_list = -1*numpy.inf * numpy.ones(matrix.shape[0], dtype=numpy.float32)
    assignments = -1 * numpy.ones(matrix.shape[0], dtype=int)
    distance_cutoff=float(distance_cutoff)

    seed=0
    generator_indices = []

    for i in xrange(cluster_cutoff):
        new_ind = seed if i == 0 else numpy.argmin(distance_list)
        print "K-centers: Finding generator %i. Will finish when % .4f is above % .4f" % (i, float(distance_list[new_ind]), float(distance_cutoff))
        if distance_list[new_ind] > distance_cutoff:
            break
        new_distance_list = matrix[new_ind, :]
        updated_indices = numpy.where(new_distance_list > distance_list)[0]
        if updated_indices.size==0:
            break
        distance_list[updated_indices] = new_distance_list[updated_indices]
        assignments[updated_indices] = new_ind
        generator_indices.append(new_ind)
    return numpy.array(generator_indices), numpy.array(assignments), numpy.array(distance_list)

                                                           
