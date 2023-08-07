library(rphast)
args = commandArgs(trailingOnly=TRUE)

chrid <- args[1]

locations <- read.csv("gene140_location.csv", header=T)

align_list <- list()
print(paste("chr",chrid,".maf", sep=""))

location_chr <- locations[locations$Chromosome.scaffold.name==chrid,]
n_loc <- nrow(location_chr)
print(location_chr)

chr <- read.msa(paste("chr",chrid,".maf", sep=""), format="MAF", seqnames=c("hg38", "gorGor5", "panTro5", "panPan2", "ponAbe2", "nomLeu3", "macNem1","calJac3"), ordered=TRUE)

for (i in 1:n_loc) {
    align_list[[i]] <- sub.msa(chr, start.col=location_chr$Gene.start..bp..10k[[i]], end.col=location_chr$Gene.end..bp..10k[[i]], refseq="hg38")
    print(location_chr$Gene.name[[i]])
    write.msa(align_list[[i]], paste(location_chr$Gene.name[[i]],"_new_8species.fasta", sep=""), format="FASTA")
}
