
### Script from adapted from
### https://elearning.bits.vib.be/courses/alphafold/ by Jasper Zuallaert (VIB-UGent), with the help of Alexander Botzki (VIB) and Kenneth Hoste (UGent).
### https://raw.githubusercontent.com/jasperzuallaert/VIBFold/main/visualize_alphafold_results.py

import glob
import math
import os
import numpy as np
from matplotlib import pyplot as plt
import sys
import pickle
from Bio.PDB import *

def generate_output_images(feature_dict, out_dir):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    # retrieving length of chains from the pdb
    clen=[]
    r=0
    parser = PDBParser()
    structure = parser.get_structure("P", pdb)
    for chain in structure[0]:
        for residue in chain.get_residues():
            r += 1
        clen.append(r)
    ##################################################################
    plt.figure(dpi=300)
    ##################################################################
    plt.title("Sequence coverage")
    plt.imshow(final,
	       interpolation='nearest', aspect='auto',
	       cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    for cl in clen: plt.axvline(cl, color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/coverage.png")

    ##################################################################

input_dir = str(sys.argv[1]) ### takes the path to features.pkl and reanked_0.pdb

feature_dict = pickle.load(open(f'{input_dir}/features.pkl','rb'))
pdb = f'{input_dir}/ranked_0.pdb'

generate_output_images(feature_dict, input_dir)
