# *POKEMON*
*pokemon is a structure based variance component test for studying rare variants*

## Setup
### Dependencies
python > 3.6  
fastlmm
  see the instruction guide here: https://pypi.org/project/fastlmm/
pandas  
  see the instruction guide here:https://pandas.pydata.org/docs/getting_started/install.html

### Installation
git clone https://github.com/bushlab-genomics/POKEMON.git 
cd POKEMON 

### Example:
gene=ENST00000373113
python run_pokemon.py --gene_name ${gene} --genotype ${gene}.raw --phenotype test.pheno --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --alpha 0 --use_aa --annotation ENST00000373113.csq --out_file results

### Flags:
--gene_name: ensemble ID, used for mapping snp from gene to protein 
--genotype: plink output with recode A option(The 'transpose' modifier).  
    The columns for genotype file is FID IID PAT MAT SEX PHENOTYPE <snp1> ... <snp2>  
    * snp must be named as chr:pos:alt:ref(e.g., 6:41129275:G:C)  
--cov_file: covariate file. *optional*  
  the columns for covariate file are: FID IID <cov1> ... <cov2> 
--cov_list: covariates to be used. *optional, but compulsory if --cov_file is used* 
  covariates name   
