#!/usr/bin/python

import scipy.signal
import scipy.stats
import numpy
import os
import argparse
import warnings
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROSETTADB = '/home/dimaio/Rosetta/main/database'
ATOMFILE = '%s/chemical/atom_type_sets/fa_standard/atom_properties.txt'%ROSETTADB
DBPATHS = [
    '%s/chemical/residue_type_sets/fa_standard/residue_types/l-caa'%ROSETTADB,
    '%s/chemical/residue_type_sets/fa_standard/residue_types/nucleic/dna'%ROSETTADB
]

SKIPATOMS = ['Hapo','Hpol','Haro','HS', 'HNbb']

BINS = numpy.arange(1.5,5.05,0.05)
P_PSEUDO = 0.01
MINCOUNT = 100


atom_dtype = numpy.dtype( [
    ("atomname", numpy.unicode_, 4),
    ("atomtype", numpy.int32),
    ("resname", numpy.unicode_, 3),
    ("chnid", numpy.unicode_, 1),
    ("resid", numpy.int32),
    ("X", numpy.float64, 3),
] )

def parse_pdb(pdbfile, aamap, allatms):
    allatoms = []
    with open(pdbfile) as pdbin:
        for line in pdbin:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if ((resname in aamap) and (atomname in aamap[resname])):
                    atmtype = allatms.index(aamap[resname][atomname])
                    split_line = (
                        #line[12:16], line[17:20], line[21], line[22:26], 
                        atomname, atmtype, resname, line[21], line[22:26], 
                        (line[30:38], line[38:46], line[46:54])
                    )
                    allatoms.append(split_line)
    return (numpy.array(allatoms, dtype=atom_dtype))

def read_params(paramfile):
    all_aas = set()
    all_atms = set()
    aamap = dict()
    with open(paramfile) as paramsin:
        aaname = ''
        for line in paramsin:
            if line.startswith('IO_STRING '):
                fields = line.rstrip().split()
                aaname = fields[1]
                all_aas.add(aaname)
                aamap[aaname] = dict()
            elif line.startswith('ATOM '):
                fields = line.rstrip().split()
                atmname = fields[1]
                atype = fields[2]
                if (atype not in SKIPATOMS):
                    all_atms.add(atype)
                    aamap[aaname][atmname] = atype

    return aamap,all_aas,all_atms


def load_db():
    all_aas = set()
    all_atms = set()
    aamap = dict()

    for dbpath in DBPATHS:
        for root, dirs, files in os.walk(dbpath):
            for file in files:
                if (file[0] == '.'):
                    continue
                resfile = os.path.join(root, file)
                aamap_i, aas_i, atoms_i = read_params(resfile)
                all_aas.update(aas_i)
                all_atms.update(atoms_i)
                aamap.update(aamap_i)
    return aamap,list(all_aas),list(all_atms)


def load_radii(all_atms):
    radii = numpy.zeros(len(all_atms))
    with open(ATOMFILE) as instream:
        for line in instream:
            words = line[:-1].split()
            atype = words[0]
            if atype in all_atms:
                radius = float(words[2])
                ai = all_atms.index(atype)
                radii[ai] = radius
    return radii

    
def get_neighbors(pdbdata):
    MAXVAL=numpy.max(BINS)
    minX = numpy.min(pdbdata['X'],axis=0)
    hashbins = (
        numpy.floor((pdbdata['X'] - minX[None,:])/(MAXVAL+1.0))
    ).astype(numpy.int64)
    FACT = 100
    indices = FACT*FACT*hashbins[:,0]+FACT*hashbins[:,1]+hashbins[:,2]

    nbins = 0
    maxbin = 0
    pthash = dict()
    for i,ind in enumerate(indices):
        if (not ind in pthash):
            nbins += 1
            pthash[ind] = numpy.array([], dtype=numpy.int32)
        pthash[ind] = numpy.concatenate( (pthash[ind], numpy.array([i])) )
        if (len(pthash[ind]) > maxbin):
            maxbin = len(pthash[ind])
        
    halfplane = [
        0, 1,
        FACT-1, FACT, FACT+1,
        FACT*FACT-FACT-1, FACT*FACT-FACT, FACT*FACT-FACT+1,
        FACT*FACT-1, FACT*FACT, FACT*FACT+1,
        FACT*FACT+FACT-1, FACT*FACT+FACT, FACT*FACT+FACT+1,
    ]

    MAXLEN = 14 * nbins * maxbin * maxbin
    dists = numpy.zeros(MAXLEN)
    types1 = numpy.zeros(MAXLEN, dtype=numpy.int32)
    types2 = numpy.zeros(MAXLEN, dtype=numpy.int32)

    outarray_i = 0
    for bin1 in pthash.keys():
        idxs1 = pthash[bin1]

        for offset in halfplane:
            bin2 = bin1+offset
            if (bin2 in pthash):
                idxs2 = pthash[bin2]
                d_ij = numpy.linalg.norm(
                    pdbdata[idxs1]['X'][:,None,:] - pdbdata[idxs2]['X'][None,:,:], axis=-1 )
                mask = numpy.logical_and(
                    d_ij<MAXVAL,
                    numpy.logical_or(
                        pdbdata[idxs1]['chnid'][:,None] != pdbdata[idxs2]['chnid'][None,:],
                        numpy.abs(pdbdata[idxs1]['resid'][:,None] - pdbdata[idxs2]['resid'][None,:]) > 9
                    )
                )
                if (offset == 0):
                    mask = numpy.tril(mask, -1)
                mask = mask.nonzero()
                nmask = len(mask[0])

                dists[outarray_i:(outarray_i+nmask)] = d_ij[mask]
                types1[outarray_i:(outarray_i+nmask)] = pdbdata[idxs1[mask[0]]]['atomtype']
                types2[outarray_i:(outarray_i+nmask)] = pdbdata[idxs2[mask[1]]]['atomtype']

                outarray_i += nmask
    return dists[:outarray_i],types1[:outarray_i],types2[:outarray_i]


def main(files, reffile, radius, genref, plot, verbose, chainids):
    aamap,all_aas,all_atms = load_db()
    radii = load_radii(all_atms)
    natms = len(all_atms)
    dlist = {}
    

    for pdbfile in files:
        pdbdata = parse_pdb(pdbfile,aamap,all_atms)
        if (len(chainids) > 0):
            pdbdata = numpy.array(
                [ai for ai in pdbdata if ai['chnid'] in chainids], 
                dtype=atom_dtype)
        ds,a1,a2 = get_neighbors(pdbdata)
        vdwsum = radii[a1]+radii[a2]+radius
        ds,a1,a2 = ds[ds<vdwsum],a1[ds<vdwsum],a2[ds<vdwsum]
        
        for x in range(len(ds)):
            ai,aj = int(a1[x]),int(a2[x])
            if (ai<aj):
                tag = str(ai)+' '+str(aj)
            else:
                tag = str(aj)+' '+str(ai)
            if (tag not in dlist):
                dlist[tag] = []
            dlist[tag].append(ds[x])

    hist = {}
    counts = {}
    for x,y in dlist.items():
        counts[x] = len(y)
        if (counts[x] > 1):
            kde = scipy.stats.gaussian_kde(y, 'scott')
            hist[x] = kde(BINS)

    if (genref):
        toplot = {}
        with open(reffile,'w') as outfile:
            for x,y in hist.items():
                cij = int(counts[x])
                i,j = [int(k) for k in x.split()]
                outfile.write("{} {} {}".format(all_atms[i],all_atms[j],cij))
                for y_i in y:
                    outfile.write(" {:.6f}".format(y_i))
                outfile.write("\n")

                if (plot):
                    name = all_atms[i]+':'+all_atms[j]+'\n'+str(cij)
                    toplot[name] = y

        if (plot):
            nplots = len(toplot.keys())
            fig, axs = plt.subplots(int((nplots+9)/10), 10, figsize=[40,40], frameon=False)
            for count,(name,ref) in enumerate(toplot.items()):
                i,j = int(count/10),count%10
                axs[i,j].set_title(name)
                axs[i,j].plot(ref,'r')
                axs[i,j].set_yticklabels([])
                axs[i,j].set_xticklabels([])
            plt.savefig('ref.pdf')  
            
    else:
        # compare to ref
        refhist, refcounts = dict(), dict()
        with open(reffile) as infile:
            for line in infile:
                fields = line.split()
                if (int(fields[2]) >= MINCOUNT):
                    refhist[ (fields[0],fields[1]) ] = numpy.array([float(x) for x in fields[3:]])
                    refcounts[ (fields[0],fields[1]) ] = int(fields[2])

        KLsum, KLwtsum = 0.0, 0.0
        toplot = {}
        for (a1,a2) in refhist.keys():
            i1 = all_atms.index(a1)
            i2 = all_atms.index(a2)
            if (i1<i2):
                tag = str(i1)+' '+str(i2)
            else:
                tag = str(i2)+' '+str(i1)
            kldiv = (BINS[1]-BINS[0]) * numpy.sum(scipy.special.rel_entr(hist[tag]+P_PSEUDO, refhist[a1,a2]+P_PSEUDO))
            if (verbose):
                print (a1,a2,kldiv)
                
            # 
            weight = numpy.sqrt( refcounts[a1,a2] )

            KLwtsum += weight
            KLsum += weight * kldiv * kldiv

            if (plot):
                name = a1+':'+a2+'('+str(refcounts[a1,a2])+')\n'+str(kldiv)
                toplot[name] = (refhist[a1,a2], hist[tag])

        if (plot):
            nplots = len(toplot.keys())
            fig, axs = plt.subplots(int((nplots+9)/10), 10, figsize=[40,40], frameon=False)
            for count,(name,(d1,d2)) in enumerate(toplot.items()):
                i,j = int(count/10),count%10
                axs[i,j].set_title(name)
                axs[i,j].plot(d1,'r')
                axs[i,j].plot(d2,'g')
                axs[i,j].set_yticklabels([])
                axs[i,j].set_xticklabels([])
            plt.savefig('distrs.pdf')  

        ndstrs = len(refhist.keys()) 
        print('%f %d'%(sqrt(KLsum/KLwtsum), ndstrs))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='collect ssm results')
    parser.add_argument('inpdb', help='outfiles', nargs='+')
    parser.add_argument('--gen',  help='Generate reference file', action='store_true')
    parser.add_argument('--ref',  help='Reference file name', default='distr_new.REF')
    parser.add_argument('--plot',  help='Plot distributions', action='store_true')
    parser.add_argument('--verbose',  help='Be verbose', action='store_true')
    parser.add_argument('--chains',  help='Chain IDs', default='')
    parser.add_argument('--radius',  help='Radius of interactions (w.r.t. vdw radii)', default=0.5)
    args = parser.parse_args()

    chains = args.chains
    if (chains != ''):
        chains = chains.split(',')
    else:
        chains = []

    main(
        args.inpdb,
        args.ref,
        args.radius,
        args.gen,
        args.plot,
        args.verbose,
        chains
    )
