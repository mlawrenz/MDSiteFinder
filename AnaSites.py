import numpy, pylab
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

def get_residue_similarity(data, ind):
    similarity=numpy.zeros(len(data))
    for site in data:
        overlap=numpy.intersect1d(data[ind], site)
        similarity[n]=len(overlap)/len(data[ind])
    return similarity_scores

def Cluster(residue_data, cutoff=0.8, Seed=0):
    """Feed in Data with indices already selected"""
    GeneratorIndices=[Seed]
    n0,n1=residue_data.shape
    List=numpy.ones(n0)*numpy.inf
    for k in xrange(100000-1):
        print("Finding Generator %d"%(k+1))
        NewList=get_residue_similarity(residue_data, GeneratorIndices[k])
        List[numpy.where(NewList<List)]=NewList[numpy.where(NewList<List)]
        NewInd=np.argmax(List)
        if List[NewInd] < cutoff:
            break #  The current generators are good enough; do not add another one.
        GeneratorIndices.append(NewInd)
    return GeneratorIndices

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
    import pdb
    pdb.set_trace()
    residue_sites=-1*numpy.ones((numframes*num_pockets, max_num_residues))
    count=0
    for frame in xrange(numframes):
        for site in range(1, num_pockets):
            key='siteuid_%s' % site
            res=[]
            if data[key][frame]!=None:
                res=[int(i) for i in data[key][frame].split()]
                residue_sites[count][0:len(res)]=numpy.array(res)
                count+=1
            else:
                count+=1
    gens=Cluster(residue_sites)

    #plot_pocket_props(data, num_pockets)


def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-i','--input',dest='input',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":
    (options,args)=parse_cmdln()
    main(input=options.input)

