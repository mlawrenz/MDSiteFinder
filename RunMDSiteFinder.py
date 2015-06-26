import Site3D
import mdtraj
import optparse
import pylab
from numpy import * 
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

#def main(modeldir, genfile, ligandfile, proteinfile, topology, writefree=False):
def main(pocketfile, trajfile, topo):
    # load traj
    dir=os.path.dirname(trajfile)
    filename=trajfile.split('%s/' % dir)[1].split('.')[0]
    traj=mdtraj.load(trajfile, top=topo)
    newcoors=10*traj.xyz
    total_frames=newcoors.shape[0]
    # get pocket spheres, map to grid
    pocket_data=Site3D.parse_pocket_file(pocketfile)
    print "getting min max"
    reduced_coors, x_range, y_range, z_range, box_volume=Site3D.get_pocket_minmax(pocket_data, newcoors)
    space=Site3D.Site3D(xaxis=x_range, yaxis=y_range, zaxis=z_range)
    #get frew
    print "getting tally"
    tally=space.map_sphere_occupancy_grid(pocket_data, reduced_coors)
    freq=tally/total_frames
    print freq.min(), freq.max()
    space.write_dx(freq, dir, filename)


def parse_commandline():
    parser = optparse.OptionParser()
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
    main(pocketfile=options.pocketfile, trajfile=options.trajfile, topo=options.topo)


