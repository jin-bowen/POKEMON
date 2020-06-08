#import pymol
#from pymol import cmd

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import sys

def plot(var,pdb_id, chain=None): 
	pymol.finish_launching()

	if chain:
		pdb_id = pdb_id + chain

	cmd.fetch(pdb_id)
	cmd.alter(pdb_id, 'b = 0.5')
	cmd.show_as('cartoon',pdb_id)
	cmd.color('white',pdb_id)

	for i, row in var.iterrows():
		resi  = row['seq']
		chain = row['chain']
		pheno = row['pheno']
		cmd.alter('resi %s and chain %s'%(resi, chain), 'b=%s'%pheno)	  	
 	
	cmd.spectrum("b", "blue_white_red", "%s"%pdb_id, maximum=1.0, minimum=0)
	cmd.zoom()

if __name__ == "__main__":
	pdb_id = sys.argv[1]
	pdb_combo_file  = sys.argv[2]

	pdb_combo = open(pdb_combo_file, 'rb')

        genotype = pickle.load(pdb_combo)
        freqs    = pickle.load(pdb_combo)
        pheno    = pickle.load(pdb_combo)
        snps2aa  = pickle.load(pdb_combo)

	## add phenotype to the variant dataframe
	for i, row in snps2aa.iterrows():
		gd_id   = row.name
		subset  = genotype[genotype[gd_id] > 0.0]
		carrier = subset.index.values
		ipheno  = np.mean(pheno[carrier].values)
		snps2aa.at[i,'pheno'] = ipheno
	
	#plot(snps2aa, pdb_id)

	pheno_vec = snps2aa['pheno'].values.reshape(-1,1)
	pheno_mat = np.dot(pheno_vec, pheno_vec.T)

	snp = range(1, freqs.shape[0] + 1 )

	fig, ax = plt.subplots()
	im = ax.imshow(pheno_mat)
	fig.colorbar(im, ax = ax)

	#We want to show all ticks...
	ax.set_xticks(np.arange(len(snp)))
	ax.set_yticks(np.arange(len(snp)))
	# ... and label them with the respective list entries
	ax.set_xticklabels(snp)
	ax.set_yticklabels(snp)

	plt.show()




