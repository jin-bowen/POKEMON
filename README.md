# POKEMON
pokemon is a structure based variance component test for studying rare variants

## Setup
**Highly recommand to use through docker image**  
dockerhub repo: docker.io/bushlab/pokemon:latest
```bash
dir=/home/myhomespace
gene=ENSG00000183134
docker run -it --rm -v ${dir}:${dir} -w=/POKEMON docker.io/bushlab/pokemon \
python run_pokemon.py --gene_name ${gene} --genotype ${gene}.raw --phenotype test.pheno --annotation ${gene}.csq \
                     --cov_file test.cov --cov_list APOE4_dose,APOE2_dose,PC1,PC2 --alpha 0.0 --use_blosum  \
                     --out_file ${dir}/results --figures --out_fig_dir=${dir}  
```

## Dependencies
python > 3.7  

gcc > 6.3.0

Required package:(can be installed via pip)
- fastlmmclib
see the instruction guide here: https://pypi.org/project/fastlmmclib/  
- biopython
see the instruction here: https://biopython.org/wiki/Download  
- sklearn
see the instruction here: https://scikit-learn.org/stable/install.html  

Optional pacakges if visualization is on:  
-  pymol   
pymol must be installed to use the flag --figures, see instruction here: https://pymol.org/2/  

```bash
git clone https://github.com/bushlab-genomics/POKEMON.git  
cd POKEMON 
```

## Example:
```bash
gene=ENSG00000183134
python run_pokemon.py --gene_name ${gene} --genotype ${gene}.raw --phenotype test.pheno --annotation ${gene}.csq \
                     --cov_file test.cov --cov_list APOE4_dose,APOE2_dose,PC1,PC2 --alpha 0.0 --use_blosum  \
                     --out_file results 
```
## Flags:
**--gene_name**: required  
   A label used as gene name in the output result file 
   
**--genotype**: required  
   allele count matrix
   **note1: snp must be named as chr:pos:alt:ref (e.g., 6:41129275:G:C)**   
   **note2: snps must be unique**  
   
   A example for genotype file is      
| FID  | IID | PAT | MAT | SEX | PHENOTYPE | 6:41129275:G:C_G | ... | ... |
| --- | --- | --- | --- |--- | --- | --- | --- | --- |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |   

   A typical command to generate the genotype file:
```bash
plink2 --vcf <vcf file with genotypes> --set-all-var-ids @:#:\$r:\$a --snps-only --mac 1 --export A --out test_gene
```

**--phenotype**: required  
   phenotype file  
   **note: phenotypes must be unique**  
   
   format1 for testing with single phenotype: column1 for individuals, column2 for phenotype, seperated by space. 
| ID  | pheno |
| --- | ---- |
| sample1  | 0  |
| sample1  | 1  | 
| ...  | ...  |
| sample1000  | 0  |  

   format2 for testing with multiple phenotypes: row for phenotypes, columns for individuals, seperated by space.   
| pheno  | sample1 | sample2 | ... | sample1000 |
| --- | --- |--- |--- |---|
| pheno1  | 1  | 0 | ... | 1 |
| pheno2  | 1  | 0 | ... | 1 |
| ...  | ...  | ... | ... | ... |
| pheno50  | 1  | 0 | ... | 1 |


 
**--annotation**: required  
    Consequence annotations from Ensembl VEP __with vcf format__  
    INFO columns must contains CANONICAL|SWISSPROT|Amino_acids|Protein_position(can be easily achieved when run vep with outputing everything)    
    A typical script to generate the annotation file:  

```    
${dir_to_vep}/vep -i <vcf file with genotype> --cache --offline --dir <dir to cache> \
--check_ref --symbol --protein --uniprot --domains --canonical --biotype --coding_only --assembly GRCh37 \
--vcf -o test_gene.csq --no_stats
```
    
**--alpha**:  required    
    alpha = 0: using structural kernel only  
    alpha = 0.5: using combined kernel of frequency and structure  
    alpha = 1: using frequency kernel only(equivalent to a standard SKAT test)  

**--out_file**: required   
    output file where the POKEMON will write  

**--cov_file**: *optional*
    covariate file.  
    | ID  | cov1 | cov2 | ... | cov1000 |
| --- | --- |--- |--- |---|
| sample1  | 1  | 0 | ... | 1 |
| sample2 | 1  | 0 | ... | 1 |
| ...  | ...  | ... | ... | ... |
| sample1000  | 1  | 0 | ... | 1 |

**--cov_list**: *optional, but required if --cov_file is flagged*   
    covariates to be used, seperate by comma    
    **covariate must be present in the columns for covariate file**  

**--pdb**: *optional*   
    e.g., --pdb 5eli, POKEMON will run on the specificed protein:5eli rather the optimal one  
   
**--maf**: *optional*, default as 0.05    
    e.g., --maf 0.01, POKEMON will only run on variants with minor allele frequency < 0.01  
  
**--database**: *optional*, default as pdb    
    can only be pdb or alphafold     
                                                                                         
**--use_aa**: *optional*  
    if explicitly flagged, the kernel will be further scaled by AA change weight from BLOSUM62 matrix  

**--figures**: *optional*    
    if explicitly flagged, POKEMON will save pymol figures  
    **pymol must be installed to use this flag**  
  
**--out_fig_dir**: *optional but required if --figures is flagged*   
    Directory where POKEMON will write figures to                                                                                         

### Reference:  
Jin, B., Capra, J.A., Benchek, P., Wheeler, N., Naj, A.C., Hamilton-Nelson, K.L., Farrell, J.J., Leung, Y.Y., Kunkle, B., Vadarajan, B., et al. (2021). An Association Test of the Spatial Distribution of Rare Missense Variants within Protein Structures Improves Statistical Power of Sequencing Studies. BioRxiv 2021.08.09.455695.    

### License:  
MIT   
