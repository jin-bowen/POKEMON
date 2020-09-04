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
python run_pokemon.py --gene_name ENST00000373178  --genotype test.vcf --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --alpha 1.0  --ref_mapping ref.tab --ref_pdb pdb_ca_coord/ --out_file test
