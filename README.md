# Overview 
pokemon is a structure based variance component test for studying rare variants

# Dependencies
python > 3.6  
fastlmm  
pandas  

# Installation
git clone https://github.com/bushlab-genomics/POKEMON.git 

cd $HOME/POKEMON 

# Example:
gene=ENST00000373113
python run_pokemon.py --gene_name ${gene} --genotype ${gene}.raw --phenotype test.pheno --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --alpha 0 --use_aa --annotation ENST00000373113.csq --out_file results

