# POKEMON
pokemon is a structure based variance component test for studying rare variants

## Setup
### Dependencies
- python > 3.6  
- fastlmm  
see the instruction guide here: https://pypi.org/project/fastlmm/

### Installation
```bash
git clone https://github.com/bushlab-genomics/POKEMON.git  
cd POKEMON 
```
### Example:
```
gene=ENST00000373113
python run_pokemon.py --gene_name ${gene} --genotype ${gene}.raw --phenotype test.pheno --cov_file test.cov --cov_list APOE4_dose,APOE2_dose --alpha 0 --use_aa --annotation ${gene}.csq --out_file results
```
### Flags:
**--gene_name**: required  
   ensemble ID, used for mapping snp from gene to protein  
   
**--genotype**: required  
   plink output with recode A option(The 'transpose' modifier).    
   The columns for genotype file is FID IID PAT MAT SEX PHENOTYPE <snp1> ... <snp2>    
   **snp must be named as chr:pos:alt:ref (e.g., 6:41129275:G:C)**
  
   A typical command to generate the genotype file:
   ```
   bcftools annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' <vcf file with genotype>
   plink --vcf <vcf file with genotype> --snps-only  --allow-no-sex --max-maf 0.05 --recode A --threads 4 --out test_gene
   # change the formart from chr:pos:ref:alt_alt -> chr:pos:ref:alt
   sed -i 's/_[A-Z]//g' test_gene.raw
   ```
   
**--cov_file**:  *optional*   
    covariate file.  
    the columns for covariate file are: FID IID <cov1> ... <cov2>
    A typical command to generate the covariate file:
    ```
    ${dir_to_plink}/plink --vcf <vcf file with genotype> --allow-no-sex **--covar <vcf file with genotype>** --prune --snps-only  --allow-no-sex --max-maf 0.05 --recode A --threads 4 --out test_gene
    ```  
   
**--cov_list**: *optional, but compulsory if --cov_file is used*   
    covariates to be used  
    **covariate must be present in the columns for covariate file**  
 
**--annotation**: required  
    Consequence annotations from Ensembl VEP __with vcf format__  
    INFO columns must contains CANONICAL|SWISSPROT|Amino_acids|Protein_position(can be easily achieved when run vep with outputing everything)    
    A typical script to generate the annotation file:
    ```
    ${dir_to_vep}/vep -i --vcf <vcf file with genotype> --format vcf --cache --offline --dir <dir to cache> --check_existing --symbol --protein --uniprot --domains --canonical --biotype --pubmed --coding_only --assembly GRCh37 --buffer_size 50000  --fork 8 --vcf -o test_gene.csq --no_stats
    ```    
**--alpha**:  required    
    alpha = 0: using structural kernel only  
    alpha = 0.5: using combined kernel of frequency and structure  
    alpha = 1: using frequency kernel only(equivalent to a standard SKAT test)  

**--use_aa**: *optional*  
    if explicitly flagged, the kernel will be further scaled by AA change weight from BLOSUM62 matrix  
  
**--out_file**: required  
    output file where the POKEMON will write to  

### Reference:  
    TBD  

### License:  
    MIT  
