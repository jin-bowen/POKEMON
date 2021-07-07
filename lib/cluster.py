from sklearn.cluster import DBSCAN
import pymol
from pymol import cmd
import numpy as np 
import pandas as pd
import seaborn as sns
from scipy import sparse
import pickle
import sys

class cluster:

	def __init__(self,genotype,snps2aa,phenotype,distance,pdb):
		self.genotype  = genotype
		self.phenotype = phenotype
		self.snps2aa   = snps2aa
		self.distance  = distance
		self.pdb       = pdb

		self.snps2aa = self.snps2aa[self.snps2aa['structure']==self.pdb]

	def cluster(self,matrix):
		clustering = DBSCAN(metric='precomputed',eps=14,min_samples=3)
		clustering.fit(matrix)
		label = clustering.labels_
		return label

	def cluster_analysis(self):
		self.score_on_var()
		case_mtx = self.distance[self.case_id,:][:,self.case_id]
		ctrl_mtx = self.distance[self.ctrl_id,:][:,self.ctrl_id]
		# case variant
		case_label = self.cluster(case_mtx)
		# ctrl variant
		ctrl_label = self.cluster(ctrl_mtx)

		case_cluster = pd.DataFrame(data=case_label,index=self.case_id, columns=['grp'])
		self.case_cluster_df = pd.merge(case_cluster,self.snps2aa,left_index=True,right_on='id')
		ctrl_cluster = pd.DataFrame(data=ctrl_label,index=self.ctrl_id, columns=['grp'])
		self.ctrl_cluster_df = pd.merge(ctrl_cluster,self.snps2aa,left_index=True,right_on='id')

	def plot_cluster(self,outfile):
		# for case
		n_case_grp = self.case_cluster_df['grp'].max() + 1
		case_color_list = sns.color_palette("OrRd", n_case_grp)
		case_color_dict = dict(zip(list(range(n_case_grp)),case_color_list))

		#pymol.finish_launching()
		cmd.fetch(self.pdb)
		cmd.show_as('cartoon',self.pdb)
		cmd.color('white',self.pdb)
		
		self.case_cluster_df_grp = self.case_cluster_df.groupby('grp')
		for igrp, case_grp in self.case_cluster_df_grp:
			if igrp == -1: continue
			color_name = 'case'+str(igrp)
			cmd.set_color(color_name, case_color_dict[igrp]) 
			for i, row in case_grp.iterrows():
				resi  = row['structure_position']
				chain = row['chain']
				selec = 'snp%s_case%s'%(i,igrp)
				selec_atom = 'snp_atom%s_case%s'%(i,igrp)
				cmd.select(selec,'name ca and resi %s and chain %s'%(resi,chain))
				cmd.create(selec_atom, selec)
				cmd.set("sphere_scale",0.8)
				cmd.show('sphere',selec_atom)
				cmd.color(color_name, selec_atom)

		# for ctrl
		n_ctrl_grp = self.ctrl_cluster_df['grp'].max() + 1
		ctrl_color_list = sns.color_palette("YlGnBu", n_ctrl_grp)
		ctrl_color_dict = dict(zip(list(range(n_ctrl_grp)),ctrl_color_list))

		self.ctrl_cluster_df_grp = self.ctrl_cluster_df.groupby('grp')
		for igrp, ctrl_grp in self.ctrl_cluster_df_grp:
			if igrp == -1: continue
			color_name = 'ctrl'+str(igrp)
			cmd.set_color(color_name, ctrl_color_dict[igrp]) 
			for i, row in ctrl_grp.iterrows():
				resi  = row['structure_position']
				chain = row['chain']
				selec = 'snp%s_ctrl%s'%(i,igrp)
				selec_atom = 'snp_atom%s_ctrl%s'%(i,igrp)
				cmd.select(selec,'name ca and resi %s and chain %s'%(resi,chain))
				cmd.create(selec_atom, selec)
				cmd.set("sphere_scale",0.8)
				cmd.show('sphere',selec_atom)
				cmd.color(color_name,selec_atom)

		cmd.bg_color("white")
		cmd.zoom()
		cmd.png('%s_cluster.png'%outfile, dpi=300)
		cmd.save('%s_cluster.pse'%outfile)

	def plot(self,outfile):

		snp_score_scale = self.snp_score
		snp_score_scale['es'] = (snp_score_scale['es'] + 1)/2
		df = pd.merge(snp_score_scale, self.snps2aa,on='id')
	
		pymol.finish_launching()
		cmd.fetch(self.pdb)
		cmd.alter(self.pdb, 'b = 0.5')
		cmd.show_as('cartoon',self.pdb)
		cmd.color('white',self.pdb)
		
		for i, row in df.iterrows():
			resi  = row['structure_position']
			chain = row['chain']
			pheno = row['es']
			cmd.select(i, 'name ca and resi %s'%(resi))
			cmd.set("sphere_scale", 0.5)
			cmd.show('sphere', i)
			cmd.alter('resi %s and chain %s'%(resi, chain), 'b=%s'%pheno)
		
		cmd.spectrum("b", "blue_white_red", "%s"%self.pdb, maximum=1.0, minimum=0.0)
		cmd.bg_color("white")
		cmd.zoom()
		cmd.png('%s.png'%outfile, dpi=300)
		cmd.save('%s.pse'%outfile)

	def score_on_var(self):
	
		geno  = self.genotype.values
		geno_sum = geno.sum(axis=0)
		pheno = self.phenotype.values.reshape((1,-1))
		pheno = (pheno - 0.5 ) * 2  

		percent = pheno.dot(geno).reshape(-1)
		percent = percent / geno_sum

		percent_df = pd.DataFrame(index=self.genotype.columns)
		percent_df['es'] = percent
	
		snp_df = pd.merge(percent_df, self.snps2aa,left_index=True, right_on='varcode')
		snp_df_sub = snp_df.loc[snp_df['structure']==self.pdb,['id','es']]
		snp_df_sub.dropna(inplace=True)

		ctrl_id = snp_df_sub.loc[snp_df_sub['es']<0, 'id']
		case_id = snp_df_sub.loc[snp_df_sub['es']>0, 'id']

		self.ctrl_id = ctrl_id
		self.case_id = case_id
		self.snp_score = snp_df_sub
def main():

	file  = sys.argv[1]
	pdb = sys.argv[2]
	prefix=file.split('.')[0]

	with open(file, 'rb') as f:
		genotype = pickle.load(f)
		phenotype = pickle.load(f)
		snps2aa = pickle.load(f)
		distance_mat = pickle.load(f)

	cls = cluster(genotype,snps2aa,phenotype,distance_mat,pdb)
	outfig = prefix
	cls.cluster_analysis()
	cls.plot(outfig)
	cls.plot_cluster(outfig)	
	
# main body
if __name__ == "__main__":
	main()

