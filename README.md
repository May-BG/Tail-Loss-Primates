# Tail-Loss-Primates
This repository contains scripts used for identifying hominoid-specific variants in genes related to tail development. 
## Overview
### 01_identify_mutations
In order to identify variants associated with the tail loss phenotype, we conducted a comparative phylogenetic study to scan for hominoid-specific variants in 140 genes (+/-10kb) related to tail development. The full gene list is attached in file gene140_location.csv (adapted from Ensembl BIOMART). Tbxt gene is used as an example for the pipeline.  
Taking Tbxt gene as an example:  
First extract the multiple sequence alignments for each gene. Whole genome alignment can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/. Eight sequences including six hominoid species ("hg38", "gorGor5", "panTro5", "panPan2", "ponAbe2", "nomLeu3") and  two non-hominoid species("macNem1","calJac3") are extracted. Example output is TBXT_new_8species.fasta
```
Rscript 140gene_fasta_new_8species.R "$chromosomeID"  
```
Then hominoid-specific variants can be identified by
```
python mutation.py TBXT_new_8species.fasta   
```
For 140 genes, after gathering fasta file for each gene, run:  
```
bash forloop_python.sh
```
### 02_classify_mutations
inputs/ contain 140 csv files generated from step 1 (identify mutations). For each of the 140 genes, variants are classified into insertions, deletions and SNPs. Each type of variant is stored as a bed file in the outputs/. 
```
Rscript mutation_classifier.R
```
### 03_predict_via_vep
First concatenate all the SNPs, deletions and insertions among 140 genes
```
cat 02_classify_mutations/*_snp.bed > total_snp.bed  
cat 02_classify_mutations/*_deletion.bed > total_del.bed  
cat 02_classify_mutations/*_insertion.bed > total_ins.bed
```
Outputs are saved in  03_predict_via_vep/.
### 04_filter_vep_results
```
Rscript filter_vep_visualization.R
```
Output files and plot are stored in 04_filter_vep_results/
