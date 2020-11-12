import pymol
from pymol import cmd
import numpy as np 
import pandas as pd

def plot(var, pdb_id, dir):

	#pymol.finish_launching()
	cmd.fetch(pdb_id)
	cmd.alter(pdb_id, 'b = 0.5')
	cmd.show_as('cartoon',pdb_id)
	cmd.color('white',pdb_id)

	for i, row in var.iterrows():
		resi  = row['structure_position']
		chain = row['chain']
		pheno = row['es'] 
		cmd.select(i, 'name ca and resi %s'%(resi))
		cmd.set("sphere_scale", 0.5)
		cmd.show('sphere', i)
		cmd.alter('resi %s and chain %s'%(resi, chain), 'b=%s'%pheno)

	cmd.spectrum("b", "blue_white_red", "%s"%pdb_id, maximum=1.0, minimum=0.0)
	cmd.bg_color("white")
	cmd.zoom()
	cmd.png('%s/%s.png'%(dir,pdb_id), dpi=300)
	cmd.save('%s/%s.pse'%(dir,pdb_id))

def score_on_var(genotype,snps2aa,phenotype,pdb,dir):

	a = genotype.sum(axis=0)
	b = pd.merge(a.to_frame(), snps2aa, left_index=True, right_on=['varcode'])

	geno  = genotype.values
	geno_sum = geno.sum(axis=0)
	pheno = phenotype.values.reshape((1,-1))
	score = pheno.dot(geno).reshape(-1)
	score = score / geno_sum
	# normalize to [0,1]
	score = score - min(score) / (max(score) - min(score))
	
	score_df = pd.DataFrame(index=genotype.columns)
	score_df['es'] = score

	snp_df = pd.merge(score_df, snps2aa,left_index=True, right_on='varcode')
	snp_df_sub = snp_df[snp_df['structure']==pdb]

	plot(snp_df_sub, pdb, dir)

# main body
if __name__ == "__main__":
	main()

