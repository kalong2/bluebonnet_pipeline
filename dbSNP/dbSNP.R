args = commandArgs(trailingOnly=TRUE)

#biomart
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#library("biomaRt", quietly = T)
#snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#nt.biomart <- getBM(c("refsnp_id","chr_name","chrom_start","chrom_end", "allele"), filters="snp_filter", values=args[1], mart=snp.db)
#lastrow=as.data.frame(tail(nt.biomart, n=1))
#names=c("rsid","chrom","nuc_start","nuc_end","nuc_allele")
#for(i in 1:ncol(lastrow)) {
#	out=paste(names[i], lastrow[,i], "\n")
#	cat(out)
#}

#rsnps
#install.packages("rsnps")
library("rsnps")
out<-ncbi_snp_query(args[1])
out<-as.data.frame(out)
myvars<-c("maf","maf_population")
filtered<-out[myvars]
as.list(filtered)
