import numpy, pylab
import ClusterTanimotoMatrix as ct
import optparse


def make_numpy_array(data, key_types):
    # convert to numpy arrays and reorder
    for key in data.keys():
        if key_types[key]!='float':
            continue
        tmp=numpy.zeros((len(data[key])), dtype=key_types[key])
        for (n, entry) in enumerate(data[key]):
            try:
                tmp[n]=float(entry)
            except TypeError:
                tmp[n]=-1000
        data[key]=tmp
    return data

def get_key_types(data):
    key_types=dict()
    num_pockets=0
    max_num_residues=0
    for key in data.keys():
        if 'uid' in key:
            type='str'
            max_data=numpy.max([len(i.split()) for i in data[key] if i!=None])
            if max_data > max_num_residues:
                max_num_residues=max_data
        if 'hyd' in key:
            type='float'
        if 'res' in key:
            type='str' 
        if 'plb' in key:
            type='float'
        if 'size' in key:
            type='float'
            num_pockets+=1
        key_types[key]=type
    return key_types, num_pockets, max_num_residues

def plot_pocket_props(data, num_pockets, reorder=False):
    for site in range(1, num_pockets):
        key='sitesize_%s' % site
        percent_open=100*float(len(numpy.where(data[key]>0)[0]))/len(data[key])
        key='siteplb_%s' % site
        plb_avg=numpy.mean(data[key])
        key='sitesize_%s' % site
        if percent_open > 20.0 and plb_avg > 0.5:
            print "site %s %% open %s" % (site, percent_open)
            pylab.figure()
            if reorder==True:
                pylab.scatter(range(0, len(frame_order)), data[key][frame_order]) #, c=colors[n])
            else:
                pylab.scatter(range(0, len(data[key])), data[key], label='open %s%%' % int(percent_open))
            pylab.ylim(0, max(data[key])+10)
            pylab.legend(scatterpoints=1)
            pylab.title(key)
            key='siteplb_%s' % site
            pylab.figure()
            if reorder==True:
                pylab.scatter(range(0, len(frame_order)), data[key][frame_order]) #, c=colors[n])
            else:
                pylab.scatter(range(0, len(data[key])), data[key]) #, c=colors[n])
            pylab.ylim(0, 4)
            pylab.title(key)
    pylab.show()

def get_residue_similarity(data1, data2):
    frames=numpy.where(data1!=-1)[0]
    if not frames.size:
        sim=-1
    else:
        overlap=numpy.intersect1d(data1[numpy.where(data1!=-1)], data2[numpy.where(data2!=-1)])
        sim=float(len(overlap))/len(data1[numpy.where(data1!=-1)])
    return sim

def get_site_similarity(data, numframes, num_pockets):
    residue_sites=[]
    count=0
    map_sites=dict()
    for frame in xrange(numframes):
        for site in range(1, num_pockets):
            key='siteuid_%s' % site
            if data[key][frame]!=None:
                res=[int(i) for i in data[key][frame].split()]
                residue_sites.append(numpy.array(res))
                map_sites[count]=(frame, site)
                count+=1
            else:
                pass
    similarity_matrix=-1*numpy.ones((len(residue_sites), len(residue_sites)))
    for i in xrange(len(residue_sites)):
        for j in xrange(len(residue_sites)):
            sim=get_residue_similarity(residue_sites[i], residue_sites[j])
            similarity_matrix[i,j]=sim
    return residue_sites, similarity_matrix, map_sites    

def main(input, reorder=False):
    #if frames are not ordered properly in mdb file (546 after 5455, etc)
    if reorder==True:
        frame_order=numpy.loadtxt('order_of_frames', dtype=str)
        numbers=numpy.array([int(i.split('frames/frame')[1].split('.pdb')[0]) for i in frame_order])
        frame_order=numpy.argsort(numbers)
    #get entries from mdb text file
    file=open(input)
    line=file.readline()
    data=dict()
    for label in line.split(','):
        label=label.rstrip().rstrip('"').split('"')[1]
        data[label]=[]
    key_order=[i.rstrip().rstrip('"').split('"')[1] for i in line.split(',')]
    print key_order
    file=open(input)
    start=0
    numframes=len(file.readlines())-1
    file.close()
    file=open(input)
    for line in file.readlines():
        column=0
        if start==0:
            start=1
            continue
        else:
            for entry in line.split(','):
                key=key_order[column]
                if len(entry.rstrip().rstrip('"'))==0:
                    data[key].append(None)
                    column+=1
                else:
                    data[key].append(entry.rstrip().rstrip('"').split('"')[1])
                    column+=1
    key_types, num_pockets, max_num_residues=get_key_types(data)
    data=make_numpy_array(data, key_types)
    numframes=10
    residue_sites, similarity_matrix, map_sites=get_site_similarity(data, numframes, num_pockets)
    import pdb
    pdb.set_trace()
    cutoff=0.8
    gens, assignments, distances=ct.cluster(similarity_matrix, distance_cutoff=cutoff, cluster_cutoff=None)
    # assignments link to residue site list, which is linked to binding site
    # number for a certain frame
    #####
    # goal is to rename the sites according to their cluster, so site3 in frame
    # 1000 may become site 1
    # is this right??
    for (n, site) in enumerate(assignments):
        (orig_frame, orig_site)=map_sites[n]
        (new_frame, new_site)=map_sites[site]
        data[key][frame]=

    for frame in xrange(numframes):
        for site in range(1, num_pockets):
            key='siteuid_%s' % site
            if data[key][frame]!=None:
    #plot_pocket_props(data, num_pockets)


def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-i','--input',dest='input',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":
    (options,args)=parse_cmdln()
    main(input=options.input)

