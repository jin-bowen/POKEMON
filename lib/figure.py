import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.cluster import DBSCAN
import pymol
from pymol import cmd
import numpy as np 
import pandas as pd
import pickle
import sys

class figure:

	def __init__(self,genotype,snps2aa,phenotype,distance,pdb):
		self.genotype  = genotype
		self.phenotype = phenotype
		self.snps2aa   = snps2aa
		self.distance  = distance
		self.pdb       = pdb

		self.snps2aa = self.snps2aa[self.snps2aa['structure']==self.pdb]
		snps_sum = self.genotype.sum(axis=0)
		snps_sum = snps_sum[snps_sum>0]
		snps2aa_subset = self.snps2aa.merge(snps_sum.to_frame(),left_on='varcode',right_index=True)
		number_pheno_val = len(np.unique(self.phenotype.values.reshape(-1)))

		if number_pheno_val < 3:
			self.pheno_class = 'bin'
		else: self.pheno_class = 'cont'

		print('snp:', len(snps_sum))
		print('mapped snp:',snps2aa_subset.shape[0])
		print('#sub:', self.phenotype.shape[1])
		print('Protein:',self.pdb)
		
		if self.pheno_class == 'bin':
			print('#case sub:', np.sum(self.phenotype.values))
		elif self.pheno_class == 'cont':
			print("mean:", np.mean(phenotype.values.reshape(-1)))
			print("quantile:",np.quantile(phenotype.values.reshape(-1), q=0.25),\
				np.quantile(phenotype.values.reshape(-1), q=0.75))

	def cluster(self,matrix):
		clustering = DBSCAN(metric='precomputed',eps=14,min_samples=3)
		clustering.fit(matrix)
		label = clustering.labels_
		return label

	def cluster_analysis(self):
		case_id, ctrl_id, snp_df_sub = self.score_on_var()	
		case_mtx = self.distance[case_id,:][:,case_id]
		ctrl_mtx = self.distance[ctrl_id,:][:,ctrl_id]
		# case variant
		case_label = self.cluster(case_mtx)
		# ctrl variant
		ctrl_label = self.cluster(ctrl_mtx)

		case_cluster = pd.DataFrame(data=case_label,index=case_id, columns=['grp'])
		case_cluster_df = pd.merge(case_cluster,self.snps2aa,left_index=True,right_on='id')
		ctrl_cluster = pd.DataFrame(data=ctrl_label,index=ctrl_id, columns=['grp'])
		ctrl_cluster_df = pd.merge(ctrl_cluster,self.snps2aa,left_index=True,right_on='id')

		case_mtx_df = pd.merge(case_cluster_df, snp_df_sub, on='id')
		ctrl_mtx_df = pd.merge(ctrl_cluster_df, snp_df_sub, on='id')

		case_mtx_df.to_csv('%s_case.csv'%self.pdb, index=False)
		ctrl_mtx_df.to_csv('%s_ctrl.csv'%self.pdb,index=False)

		return case_cluster_df, ctrl_cluster_df
	
	def plot_cluster(self,outfile):
		case_cluster_df, ctrl_cluster_df = self.cluster_analysis()
		print(case_cluster_df)
		case_cluster_df.to_csv('%s.csv'%outfile,mode='w')
		print(ctrl_cluster_df)
		ctrl_cluster_df.to_csv('%s.csv'%outfile,mode='a')
		# for case
		n_case_grp = case_cluster_df['grp'].max() + 1
		case_map = plt.get_cmap("OrRd_r", n_case_grp + 3)
		case_color_list = [case_map(i)[:3] for i in range(n_case_grp)]
		case_color_dict = dict(zip(list(range(n_case_grp)),case_color_list))

		#pymol.finish_launching()
		cmd.reinitialize()
		try:
			cmd.load('ref/pdb/pdb%s.ent'%self.pdb, self.pdb)
		except:
			cmd.load('ref/alphafold/%s.pdb'%self.pdb, self.pdb)
		
		cmd.show_as('cartoon',self.pdb)
		cmd.color('white',self.pdb)
		
		case_cluster_df_grp = case_cluster_df.groupby('grp')
		for igrp, case_grp in case_cluster_df_grp:
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
		n_ctrl_grp = ctrl_cluster_df['grp'].max() + 1
		ctrl_map = plt.get_cmap("PuBu_r", n_ctrl_grp + 3)
		ctrl_color_list = [ctrl_map(i)[:3] for i in range(n_ctrl_grp)]
		ctrl_color_dict = dict(zip(list(range(n_ctrl_grp)),ctrl_color_list))

		ctrl_cluster_df_grp = ctrl_cluster_df.groupby('grp')
		for igrp, ctrl_grp in ctrl_cluster_df_grp:
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
				cmd.set("sphere_scale",0.6)
				cmd.show('sphere',selec_atom)
				cmd.color(color_name,selec_atom)

		cmd.bg_color("white")
		cmd.zoom()
		cmd.orient()
		cmd.png('%s_cluster.png'%outfile, width=2400, height=2400, dpi=300, ray=1)
		cmd.save('%s_cluster.pse'%outfile)

	def plot(self,outfile):

		if self.pheno_class == 'bin':	
			ctrl_id, case_id, snp_df_sub = self.score_on_var()	
		else: snp_df_sub = self.score_on_var()
		df = pd.merge(snp_df_sub, self.snps2aa,on='id')
	
		#pymol.finish_launching()
		cmd.reinitialize()
		try:
			cmd.load('ref/pdb/pdb%s.ent'%self.pdb, self.pdb)
		except:
			cmd.load('ref/alphafold/%s.pdb'%self.pdb, self.pdb)
		
		cmd.show_as('cartoon',self.pdb)
		cmd.color('white',self.pdb)

		cmap = mcolors.LinearSegmentedColormap.from_list("", ["tab:blue","white","tab:red"])
	
		if self.pheno_class == 'bin':
			df['es_norm'] = df['es']
		elif self.pheno_class == 'cont':
			df['es_norm'] = (df['es']-df['es'].min())/(df['es'].max()-df['es'].min()) 

		norm = mcolors.Normalize(vmin=0, vmax=1)
		for i, row in df.iterrows():
			resi  = row['structure_position']
			chain = row['chain']
			pheno = row['es_norm']
			selec = 'snp%s'%i
			selec_atom = 'snp_atom%s'%i

			icolor = cmap(norm(pheno))[0:3]
			color_name = 'color%s'%i
			cmd.set_color(color_name, icolor)

			cmd.select(selec,'name ca and resi %s and chain %s'%(resi,chain))
			cmd.create(selec_atom, selec)
			cmd.set("sphere_scale", 0.6)
			cmd.show('sphere', selec_atom)
			cmd.color(color_name, selec_atom)

		cmd.bg_color("white")
		cmd.zoom()
		cmd.orient()
		cmd.save('%s.pse'%outfile)
		cmd.png('%s.png'%outfile, width=2400, height=2400, dpi=300, ray=1)

	def score_on_var(self):
		geno  = self.genotype.values
		geno_sum = geno.sum(axis=0)
		pheno = self.phenotype.values.reshape((1,-1))

		percent = pheno.dot(geno).reshape(-1)
		percent = percent / geno_sum

		percent_df = pd.DataFrame(index=self.genotype.columns)
		percent_df['es'] = percent
	
		snp_df = pd.merge(percent_df, self.snps2aa,left_index=True, right_on='varcode')
		snp_df_sub = snp_df.loc[snp_df['structure']==self.pdb,['id','es']]
		snp_df_sub.dropna(inplace=True)

		if self.pheno_class == 'bin':
			ctrl_id = snp_df_sub.loc[snp_df_sub['es']<=0.5, 'id']
			case_id = snp_df_sub.loc[snp_df_sub['es']>0.5, 'id']
			return case_id, ctrl_id, snp_df_sub
		else: return snp_df_sub

# main body
if __name__ == "__main__":
	main()

