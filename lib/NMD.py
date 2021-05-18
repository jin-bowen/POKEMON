from prody import *
import numpy as np 
import pandas as pd
from modeller import *
import sys
import pickle

def mutate_structure(code,var):
	env = environ()
	# Read the topology library with non-hydrogen atoms only:
	env.libs.topology.read(file='$(LIB)/top_heav.lib')
	
	# Read the CHARMM parameter library:
	env.libs.parameters.read(file='$(LIB)/par.lib')
	# Read the original PDB file and copy its sequence to the alignment array:
	aln = alignment(env)
	mdl = model(env, file=code)
	aln.append_model(mdl, atom_files=code, align_codes=code)
	# Select the residues to be mutated: in this case all ASP residues:

	for record in var.iterrows():
		ires = str(record['structure_position'])
		target = record['alt_aa_full']
		sel = selection(mdl.residues[ires])
		sel.mutate(residue_type=target)

	# sequence in the alignment):
	mutated_code = code + '-1'
	aln.append_model(mdl, align_codes=mutated_code)
	# Generate molecular topology for the mutant:
	mdl.clear_topology()
	mdl.generate_topology(aln[mutated_code])
	# Transfer all the coordinates you can from the template native structure
	# to the mutant (this works even if the order of atoms in the native PDB
	# file is not standard):
	mdl.transfer_xyz(aln)
	# Build the remaining unknown coordinates for the mutant:
	mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
	# Write the mutant to a file:
	mdl.write(file='%s.atm'%mutated_code)


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

def main():
	
	file  = sys.argv[1]
	pdb = sys.argv[2]
	
	with open('aa.chart') as f:
		aa_chart = dict(x.rstrip().split(None, 1) for x in f)

	df = pd.read_csv(file,header=0, index_col=0)
	df[['ref_aa','alt_aa']] = df['aa'].str.split('/', expand=True)
	df['alt_aa_full'] = df['alt_aa'].apply(lambda x: aa_chart[x])
	
	code='2esb'
	mutate_structure(code,df)
	
# main body
if __name__ == "__main__":
	main()

