library(dplyr)
library(stringr)
library(tidyr)
library(R.utils)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)

location_new <- read.csv("../01_identify_mutations/gene140_location.csv", header = T)

name.list <- location_new$Gene.name
### rewrite codes to extract the sequences along with coordinate:
classify_mutation <- function(name){
  
  mutation_file_name <- paste("inputs/",name,"_8species_03.csv", sep="")
  lin28a_mutation <- read.csv(mutation_file_name, header=F)
  lin28a_mutation <- lin28a_mutation %>% filter(as.character(V2)!=as.character(V5) | as.character(V2)!=as.character(V4))

  lin28a_mutation$coordinate <- lin28a_mutation$V1+location_new[location_new$Gene.name==toupper(as.character(name)),"Gene.start..bp."]-10000
  lin28a_mutation <- lin28a_mutation[order(lin28a_mutation$V7),]

  lin28a_locations <- location_new %>% filter(as.character(Gene.name)==toupper(name)) %>% dplyr::select(Gene.start..bp., Gene.end..bp.)
  lin28a_locations$region <- name
  colnames(lin28a_locations) <- c("start","end","region")
  
    lin28a_mutation_REGIONS <- lin28a_mutation %>% 
    crossing(lin28a_locations) %>% 
    filter(coordinate >= start-10000 & coordinate <= end+10000)

  lin28a_mutation_REGIONS[order(lin28a_mutation_REGIONS$V6),] -> lin28a_mutation_REGIONS
  
  lin28a_mutation_REGIONS %>% filter(V2 %in% c("*","-","N")) ->  lin28a_mutation_REGIONS_gaps_in_human
  
  dat <- lin28a_mutation_REGIONS_gaps_in_human$V7
  
  human_deletion <- cbind(as.data.frame(seqToIntervals(dat)), seqToIntervals(dat)[,2]-seqToIntervals(dat)[,1]+1)
  colnames(human_deletion) <- c("deletion_from","deletion_to","deletion_length")
  
  del_seq <- character()
  ref <- character()
  for (i in 1:nrow(human_deletion)) {
    data <- lin28a_mutation_REGIONS_gaps_in_human %>% filter(V7 >= human_deletion[i,"deletion_from"] & V7 <= human_deletion[i,"deletion_to"]) %>% dplyr::select(V4)
    del_seq[[i]] <- paste(data$V4, collapse = "")
    ref[[i]] <- paste(rep("-", nrow(data)),collapse = "")
  }
  
  human_deletion <- cbind(human_deletion, as.data.frame(del_seq),as.data.frame(ref))

  ### combine length with data frame
  lin28a_mutation_REGIONS_gaps_in_human <- lin28a_mutation_REGIONS_gaps_in_human %>% 
    crossing(human_deletion) %>% 
    filter(V7 >= deletion_from & V7 <= deletion_to)
  
  human_deletion_wide <- merge(human_deletion, lin28a_mutation_REGIONS_gaps_in_human,  by.y = "V7", by.x = "deletion_from")
  human_deletion_trimmed <- human_deletion_wide %>% dplyr::select(region, coordinate, deletion_length.y, del_seq.y,ref.y)
  
  
  lin28a_mutation_REGIONS %>% filter(V5 %in% c("*","-","N") & V4 %in% c("*","-","N")) -> lin28a_mutation_REGIONS_insertion_in_human
  
  
  dat <- lin28a_mutation_REGIONS_insertion_in_human$coordinate
  
  
  human_insertion <- cbind(as.data.frame(seqToIntervals(dat)), seqToIntervals(dat)[,2]-seqToIntervals(dat)[,1]+1)
  colnames(human_insertion) <- c("insertion_from","insertion_to","insertion_length")
  ins_seq <- character()
  alt <- character()
  if (nrow(human_insertion)==0) {
    human_insertion_trimmed <- human_insertion  
  }
  else {
    for (i in 1:nrow(human_insertion)) {
      
      data <- lin28a_mutation_REGIONS_insertion_in_human %>% filter(coordinate >= human_insertion[i,"insertion_from"] & coordinate <= human_insertion[i,"insertion_to"]) %>% dplyr::select(V2)
      
      print(data)
      ins_seq[[i]] <- paste(data$V2, collapse = "")
      
      print(ins_seq[[i]])
      alt[[i]] <- paste(rep("-", nrow(data)), collapse = "")
    }
    
    human_insertion <- cbind(human_insertion, as.data.frame(ins_seq),as.data.frame(alt))
    human_insertion_wide <- merge(lin28a_mutation_REGIONS_insertion_in_human, human_insertion, by.x = "coordinate", by.y = "insertion_from")
    human_insertion_trimmed <- human_insertion_wide %>% dplyr::select(region, coordinate, insertion_to, insertion_length,ins_seq, alt)
  }
  
  snp <- lin28a_mutation_REGIONS %>% filter(!V5 %in% c("*","-","N") | !V4 %in% c("*","-","N")) %>% filter(!V2 %in% c("*","-","N"))
  ## write bed files for insertions
  insertion <- human_insertion_trimmed
  if (nrow(insertion) >1) {
    in_bed <- insertion %>% dplyr::select(coordinate, insertion_to,ins_seq, alt)
    in_bed$chr <- rep(as.character(location_new[location_new$Gene.name==toupper(name),"Chromosome.scaffold.name"]),nrow(in_bed))
    in_bed$type <- rep("DEL",nrow(in_bed))
    if (location_new[location_new$Gene.name==toupper(name),"Strand"]=="1") {
      strand <- "+"
    } else {
      strand <- "-"
    }
    #strand <- as.character(location_new[location_new$Gene.name==toupper(name),"Strand"])
    in_bed$strand <- rep(strand,nrow(in_bed))
    in_bed$detail <- paste(in_bed$ins_seq, in_bed$alt, sep="/")
    in_bed <- in_bed %>% dplyr::select(chr, coordinate, insertion_to, detail,strand)
    
    in_bedfile <- paste("outputs/", name, "_insertion.bed", sep="")
    write.table(in_bed, in_bedfile, sep="\t", row.names = F, quote = F, col.names =F)
  } 
  
  ## write bed files for deletions
  deletion <- human_deletion_trimmed 
  deletion$coordinate2 <- deletion$coordinate+1
  del_bed <- deletion %>% dplyr::select(coordinate, coordinate2, deletion_length.y, ref.y, del_seq.y)
  
  del_bed$chr <- rep(as.character(location_new[location_new$Gene.name==toupper(name),"Chromosome.scaffold.name"]),nrow(del_bed))
  
  
  if (location_new[location_new$Gene.name==toupper(name),"Strand"]=="1") {
    strand <- "+"
  } else {
    strand <- "-"
  }
  del_bed$strand <- rep(strand,nrow(del_bed))
  del_bed$detail <- paste(del_bed$ref.y, del_bed$del_seq.y, sep="/")
  del_bed <- del_bed %>% dplyr::select(chr, coordinate2, coordinate, detail,strand)
  
  del_bedfile <- paste("outputs/", name, "_deletion.bed", sep="")
  write.table(del_bed, del_bedfile, sep="\t", row.names = F, quote = F, col.names =F)
  
  
  ## write bed files for snps
  if (nrow(snp)>0) {
    
    
    snp_bed <- snp %>% dplyr::select(coordinate, V2, V4)
    print(head(snp_bed))
    chr <- as.character(location_new[location_new$Gene.name==toupper(name),"Chromosome.scaffold.name"])
    #print(location_new[location_new$Gene.name==toupper(name),])
    print(chr)
    snp_bed$chr <- rep(chr,nrow(snp_bed))
    if (location_new[location_new$Gene.name==toupper(name),"Strand"]=="1") {
      strand <- "+"
    } else {
      strand <- "-"
    }
    #strand <- as.character(location_new[location_new$Gene.name==toupper(name),"Strand"])
    snp_bed$strand <- rep(strand,nrow(snp_bed))
    snp_bed$end <- snp_bed$coordinate
    snp_bed$change <- paste(snp_bed$V2, "/",snp_bed$V4,sep="")
    snp_bed <- snp_bed %>% dplyr::select(chr, coordinate, end,change, strand)
    print(head(snp_bed))
    snp_bedfile <- paste("outputs/", name, "_snp.bed", sep="")
    write.table(snp_bed, snp_bedfile, sep="\t", row.names = F, quote = F, col.names =F)
  }
  
}


lapply(location_new$Gene.name, classify_mutation)

