library(dplyr)
library(stringr)
library(tidyr)
library(R.utils)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)

## read biomart outputs
biomart <- read.csv("mart_export-20.txt", sep=",")
## identify one transcript per gene based on exon number and transcript length
summ_table <- biomart %>%
  group_by(Gene.name, Transcript.stable.ID) %>%
  summarize(first_exon_number=min(Exon.rank.in.transcript),exon_number=max(Exon.rank.in.transcript), exon_row_number=n(), length=min(Transcript.length..including.UTRs.and.CDS.))
one_transcript_pergene <- summ_table %>% group_by(Gene.name) %>%
  filter(exon_row_number==max(exon_number)) %>%
  filter(length==max(length))
one_transcript_biomart <- biomart %>% filter(Transcript.stable.ID %in% one_transcript_pergene$Transcript.stable.ID)
write.csv(one_transcript_biomart, "one_transcript_biomart_new.csv",quote = F)

location_new <- read.csv("gene140_location.csv", header = T)

## read deletion vep outputs
deletion_df = read.table('deletion_vep_140genes_0711.txt',sep='\t',header=TRUE)

## filter outputs by gene names and transcriptIDs
deletion_vep <- deletion_df[deletion_df$SYMBOL %in% location_new$Gene.name,]
one_transcript_deletion_vep <- deletion_vep[deletion_vep$Feature %in% one_transcript_biomart$Transcript.stable.ID.version,]
write.csv(one_transcript_deletion_vep, 'one_transcript_deletion_vep_in_140genes_0711.csv')          
# identify deletions in coding sequences
special_deletion <- one_transcript_deletion_vep[one_transcript_deletion_vep$Consequence=="coding_sequence_variant",]
write.csv(special_deletion, "coding_deletion_0711.csv")

categories <- c("intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant", "upstream_gene_variant", "downstream_gene_variant")
one_transcript_deletion_vep$variants <- sapply(one_transcript_deletion_vep$Consequence, function(x) str_extract(x, paste(categories, collapse = "|")))

one_transcript_deletion_vep_NEW <- one_transcript_deletion_vep %>%
  mutate(variants = case_when(str_detect(variants, "coding_sequence_variant") ~ "missense_variant",
                              !str_detect(variants, "missense_variant") ~ variants
  ))
## consolidate consequences of insertions
deletion_pie_table <- as.data.frame(table(one_transcript_deletion_vep_NEW$variants))
write.csv(deletion_pie_table, "deletion_pie_140genes.csv")


## read insertion vep outputs
df = read.table('ins_vep_140genes_0711.txt',sep='\t',header=TRUE)
ins_vep <- df[df$SYMBOL %in% location_new$Gene.name,]
## filter outputs by gene names and transcriptIDs
one_transcript_ins_vep <- ins_vep[ins_vep$Feature %in% one_transcript_biomart$Transcript.stable.ID.version,]
write.csv(one_transcript_ins_vep, 'one_transcript_ins_vep_in_140genes_0711.csv')          
## visualization
categories <- c("intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "inframe_deletion","frameshift_variant", "missense_variant","coding_sequence_variant","protein_altering_variant","non_coding_transcript_exon_variant", "upstream_gene_variant", "downstream_gene_variant")
one_transcript_ins_vep$variants <- sapply(one_transcript_ins_vep$Consequence, function(x) str_extract(x, paste(categories, collapse = "|")))
one_transcript_ins_vep_NEW <- one_transcript_ins_vep %>%
  mutate(variants = case_when(str_detect(variants, "protein_altering_variant") ~ "missense_variant",
                              str_detect(variants, "missense_variant") ~ "missense_variant",
                              str_detect(variants, "frameshift_variant") ~ "missense_variant",
                              str_detect(variants, "non_coding_transcript_exon_variant") ~"other",
                              !str_detect(variants, "missense_variant") ~ variants
  ))

## filter protein altering variants
special_insertion <- one_transcript_ins_vep %>% filter(Consequence %in% c("protein_altering_variant","protein_altering_variant,NMD_transcript_variant","frameshift_variant,NMD_transcript_variant","missense_variant,coding_sequence_variant","missense_variant"))
write.csv(special_insertion, "coding_insertion_0711.csv")
## consolidate consequences of insertions
pie_table <- as.data.frame(table(one_transcript_ins_vep_NEW$variants))
write.csv(pie_table, "ins_pie_140genes.csv")

## read snp vep outputs
df = read.table('snp_vep_140gene.txt',sep='\t',header=TRUE)
snp_vep <- df[df$SYMBOL %in% location_new$Gene.name,]
## filter by transcriptID
one_transcript_snp_vep <- snp_vep[snp_vep$Feature %in% one_transcript_biomart$Transcript.stable.ID.version,]
write.csv(one_transcript_snp_vep, 'one_transcript_snp_vep_in_140genes.csv')          

## identify deleterious variants by sift and polyphen score
sift <- one_transcript_snp_vep[one_transcript_snp_vep$SIFT!="-",]
sift_del <- sift[str_detect(sift$SIFT,"deleterious"),]
sift_del_low <- sift[str_detect(sift$SIFT,"deleterious_low_confidence"),]
sift_del_high <- sift_del[!sift_del$SIFT %in% sift_del_low$SIFT,]

polyphen <- one_transcript_snp_vep[one_transcript_snp_vep$PolyPhen!="-",]
damaging <- polyphen[str_detect(polyphen$PolyPhen,"damaging"),]
PRO_damaging <- polyphen[str_detect(polyphen$PolyPhen,"probably_damaging"),]

union_relaxed <- rbind(damaging, sift_del)
union_strict <- rbind(PRO_damaging, sift_del_high)

union_relaxed[!duplicated(union_relaxed$Location),] -> union_relaxed_nonrep
union_strict[!duplicated(union_strict$Location),] -> union_strict_nonrep

write.csv(union_strict_nonrep, 'snp_union_strict_nonrep_in_140genes.csv')
write.csv(union_relaxed_nonrep, 'snp_union_relaxed_nonrep_in_140genes.csv')

## consolidate consequences of snps
snp_categories <- c("intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "missense_variant", "splice_acceptor_variant", "non_coding_transcript_exon_variant", "upstream_gene_variant", "downstream_gene_variant","splice_donor_variant","stop_gained"," stop_lost","stop_retained_variant","synonymous_variant")

one_transcript_snp_vep$variants <- sapply(one_transcript_snp_vep$Consequence, function(x) str_extract(x, paste(snp_categories, collapse = "|")))

one_transcript_snp_vep_NEW <- one_transcript_snp_vep %>%
  mutate(variants = case_when(str_detect(variants, "splice") ~ "intron_variant",
                              str_detect(variants, "stop_gained") ~ "missense_variant",
                              str_detect(variants, "stop_retained_variant") ~ "synonymous_variant",
                              str_detect(variants, "non_coding_transcript_exon_variant") ~ "other",
                              !str_detect(variants, "splice") & !str_detect(variants, "stop") &!str_detect(variants, "non_coding_transcript_exon_variant")~ variants
  ))
pie_table_snp <- as.data.frame(table(one_transcript_snp_vep_NEW$variants))
write.csv(pie_table_snp, "snp_pie_140genes.csv")

## plot
pie_table_snp_deletion <- merge(pie_table_snp, deletion_pie_table, by.x="Var1", by.y = "Var1",all = TRUE)
pie_table_snp_deletion_insertion <- merge(pie_table_snp_deletion,pie_table, by.x="Var1", by.y = "Var1",all = TRUE)
colnames(pie_table_snp_deletion_insertion) <- c("category","SNP","deletion","insertion")
pie_table_snp_deletion_insertion[is.na(pie_table_snp_deletion_insertion)] <- 0
pie_table_snp_deletion_insertion_long <- gather(pie_table_snp_deletion_insertion, variant, Counts, SNP:insertion, factor_key=TRUE)

pie_table_snp_deletion_insertion_long <- data.table(pie_table_snp_deletion_insertion_long)
pie_table_snp_deletion_insertion_long<-pie_table_snp_deletion_insertion_long[order(variant,Counts)]


## Order from large to small
pie_table_snp_deletion_insertion_long$category <-   factor(pie_table_snp_deletion_insertion_long$category,levels = pie_table_snp_deletion_insertion_long[variant=="SNP"][order(Counts)][,category])

ggplot(data=pie_table_snp_deletion_insertion_long, aes(x=Counts, y=category, fill=category)) +
  geom_text(aes(label=Counts), hjust=-0.1, color="black") +
  theme_grey(base_size=14)+
  facet_wrap(variant~.,scales = "free_x") +
  geom_bar(width = 0.5, stat = "identity") +
  theme(aspect.ratio = 2, legend.position="top") +
  scale_fill_brewer(palette="Set2") +
  labs(x="\nCounts",y=NULL,fill=NULL)

ggsave(filename = "Counts.pdf",width = 10 ,height = 8)


