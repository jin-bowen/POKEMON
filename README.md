# Overview 
pokemon is a structure based variance component test for studying rare variants

# Dependencies
python > 3.6  
fastlmm  
pandas  

# Installation
git clone https://github.com/jin-bowen/POKEMON.git  
cd $HOME/POKEMON 

# Example:
python run_pokemon.py --gene_name ADPRHL2 --genotype test.vcf --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --reference pdb_ultimate.tab
