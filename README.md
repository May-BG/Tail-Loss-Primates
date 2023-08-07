# Tail-Loss-Primates
This repository contains scripts used for identifying hominoid-specific variants in genes related to tail development. 
## Overview
### 01_identify_mutations
In order to identify variants associated with the tail loss phenotype, we conducted a comparative phylogenetic study to scan for hominoid-specific variants in 140 genes related to tail development. The full gene list is attached. Tbxt gene is used as an example for the pipeline.  
#### Taking Tbxt gene as an example:  
```
Rscript 140gene_fasta_new_8species.R "6"  
```
output includes TBXT_new_8species.fasta  
```
python mutation.py TBXT_new_8species.fasta   
```
#### For 140 genes, run:  
```
bash forloop_python.sh
```
### 02_classify_mutations
The inputs are 140 csv files generated from step01. For each of the 140 genes, variants are classified into insertion, deletion and snp. Each type is stored in a bed file. 
```
Rscript mutation_classifier.R
```
### 03_predict_via_vep
```
cat 02_classify_mutations/*_snp.bed > total_snp.bed  
cat 02_classify_mutations/*_deletion.bed > total_del.bed  
cat 02_classify_mutations/*_insertion.bed > total_ins.bed
```
### 04_filter_vep_results
