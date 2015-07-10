import Site3D
import time
import mdtraj
import optparse
import pylab
import numpy
import glob
import os
import sys
import pickle
from scipy import interpolate, integrate

"""
MDSiteFinder
===========

** To be used on MOE SiteFinder results  **

Program requires as input (see options with -h flag):
*
*

Program output is a OpenDX file SiteFrequency.dx (Occupancy for passed in site
at default 0.5 Angstrom gridspace.
This file can be loaded into VMD:
* you can run the command: vmd -pdb reference.pdb -dx SiteFrequency.dx
* make a new representation in Isosurface, and scale the isovalues ranging from
* % opened, from 0 to 1.
If your 3D surface is not centered on your protein, check that your trajectory
is aligned to the reference structure.

"""

def main(pocketfile, trajfile, topo, outname, writedx=True):
    # load traj
    resolution=0.5
    pad=7.0
    dir=os.path.dirname(trajfile)
    print "REMOVING HYDROGENS FROM CALC (REFLECTS CUTOFF)"
    # check if you already ran the calc, load framelog if so
    traj=mdtraj.load(trajfile, top=topo)
    indices=[]
    for i in traj.topology.atoms:
        if 'H' not in i.name:
            indices.append(i.index)
    newcoors=10*traj.xyz[:,indices,:] #multiply by 10 bc mdtraj scalers
    total_frames=newcoors.shape[0]
    # get pocket spheres, map to grid
    pocket_data=Site3D.parse_pocket_file(pocketfile, resolution)
    print "getting min max"
    reduced_coors, x_range, y_range, z_range, box_volume=Site3D.get_pocket_minmax(pocket_data, newcoors, pad=pad, resolution=resolution)
    space=Site3D.Site3D(total_frames, resolution=resolution, xaxis=x_range, yaxis=y_range, zaxis=z_range, reduced_coors=reduced_coors)
    #get freq
    print "getting tally"
    start=float(time.time())
    framelog=space.map_sphere_occupancy_grid(pocket_data, cutoff=3.0)
    end=float(time.time())
    elapse=end-start
    print "tallied all frames and gridpoint %0.4f sec" % elapse
    start=float(time.time())
    numpy.savetxt('%s/%s_framelog_matrix.dat' % (dir, outname), framelog)
    end=float(time.time())
    elapse=end-start
    print "saved frametally %0.4f sec" % elapse
    # work with gridpoints directly if only writing PDBs
    freq=space.pocketoccup/total_frames
    freq=numpy.round(freq, decimals=1) 
    for f in numpy.arange(0, 1.1, 0.1):
        frames=numpy.where(freq==f)[0]
        if frames.size:
            space.write_pdb(dir, outname, frames, f)
    if writedx==True:
        space.write_dx(dir, outname)


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-o', '--outname', dest='outname',
                      help='output dx name')
    parser.add_option('-y', '--topo', dest='topo',
                      help='topology')
    parser.add_option('-t', '--trajfile', dest='trajfile',
                      help='MD traj file')
    parser.add_option('-p', '--pocketfile', dest='pocketfile',
                      help='protein pocket file output from MOE SiteFinder')
    #parser.add_option('-f', action="store_true", dest="writefree", help='perform free energy calc with P-L COM distances')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    #if options.writefree==True:
    main(pocketfile=options.pocketfile, trajfile=options.trajfile, topo=options.topo, outname=options.outname)


