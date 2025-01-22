import numpy as np
import sys
import os
from Bio.PDB import *
import matplotlib.pyplot as plt
import pickle

def get_pae_plddt(model_names):
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
        basename = os.path.basename(name)
        basename = basename[basename.index('model'):]
        out[f'{basename}'] = {'plddt': d['plddt'], 'pae':d['predicted_aligned_error']}
    return out

filespath = str(sys.argv[1]) ## takes the path to .pkl files and reanked_0.pdb
print(filespath)
files = [f for f in os.listdir(filespath)]
print(files)
pklfiles =  [f for f in files if '.pkl' in f and 'result' in f]
print(pklfiles)
print('Reading list of files and monomers.')
af2pdbfile = filespath + '/ranked_0.pdb'

# reading the pdbfile
parser = PDBParser()
structure = parser.get_structure("P", af2pdbfile)
uchains = []
clen = []
for chain in structure[0]:
    uchains.append(chain.get_id())
    clen.append(len(chain))
rchainid = []
for r in range(1,1+sum(clen)):
    for c in range(1,1+len(uchains)):
        if r>sum(clen[:c-1]) and r<=sum(clen[:c]):
            rchainid.append(uchains[c-1])
# building maps
for f in range(len(pklfiles)):
    print(pklfiles[f])
    pkl = pickle.load(open(filespath+pklfiles[f],'rb'))
    pae = pkl['predicted_aligned_error']

    for c in range(len(clen)):
        y_position=sum(clen[:c])
        x_position=y_position
        plt.hlines(y_position, xmin=0, xmax=sum(clen), colors='black')
        plt.axvline(x_position, color='black')
    plt.imshow(pae, cmap="Greens_r", vmin=0, vmax=30)
    plt.colorbar(label="Expected position error (Ångströms)")
    plt.title(pklfiles[f])
    plt.xlabel("Scored residue",fontsize=16)
    plt.ylabel("Aligned residue",fontsize=16)
    figname = filespath+pklfiles[f][:-4] + '_pae.png'
    plt.savefig(figname,dpi=300)
    plt.close()
