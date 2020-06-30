#!/bin/sh

#SBATCH -J pokemon
#SBATCH -n 8 # one core
#SBATCH -N 8 # on one node
#SBATCH --mem 32000
#SBATCH --time=5-05:00:00

cat gene_list/gene${SLURM_ARRAY_TASK_ID} | while read gene
do
  chr=`cat ADSP_gene_missense/${gene}.vcf.raw | head -1 | xargs -n1 | tail -1 | awk -F":" '{print $1}' | sort | uniq`
  python run_pokemon.py --gene_name ${gene} --genotype ADSP_gene_missense/${gene}.vcf.raw --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --ref_mapping ~/pdb_ultimate/pdb_ultimate_chr${chr}.variants --ref_pdb ~/pdb_ultimate/pdb_ca_coord.txt
