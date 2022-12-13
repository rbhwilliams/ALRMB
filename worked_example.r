#--this file shows how to use ALRMB on the UASB bioreactor community data used in the manuscript
#--this file by Rohan Williams on 13 December 2022--

library(Biostrings)
library(ape)
library(kmer)
library(analogue)
library(tsne)
library(dbscan)
library(Ckmeans.1d.dp)
library(RKXM)
source(file="alrmb-functions-20221103.r")

#--read in contigs and process up to and including the formation of primary clusters
a<-readDNAStringSet(filepath="assembly_new.fasta")
a.metadata<-get.contig.metadata(a)
a.k4<-generate.kmer.frequency.matrix(a,4)
a.k4.1<-convert.kmer.frequency.to.composition(a.k4)
a.k4.1.min100kbp<-a.k4.1[rownames(a.metadata)[which(a.metadata$len>=1e5)],]
a.k4.1.min100kbp.gcd<-get.great.circle.dmat(a.k4.1.min100kbp)
a.k4.1.min100kbp.gcd.tsne<-make.tsne(a.k4.1.min100kbp.gcd)
a.k4.1.min100kbp.gcd.tsne.pc<-make.primary.clusters(a.k4.1.min100kbp.gcd.tsne,2,2,a.metadata)

#--import the coverage data
#--coverage data is here (/m/ean /p/er /b/ase /c/overage)
a.mpbc<-read.table(file="meanperbasecov_flye.lr",stringsAsFactors=F,sep="\t",header=T)
a.mpbcVec<-as.numeric(a.mpbc[,2])
names(a.mpbcVec)<-a.mpbc[,1]
#--and add to metadata
a.metadata.1<-add.coverage(a.k4.1.min100kbp.gcd.tsne.pc,a.mpbcVec)
#--so now we are ready for the computation of secondary clusters
#--please note that this /excludes/ singletons (primary clusters tagged as '0' by DBSCAN, and which may be complete genomes in their own right)
a.metadata.2<-make.secondary.clusters(a.metadata.1,"cov")

#--finally, export the contigs in each secondary cluster to FASTA for analysis using QUAST, Check and GTDB
#--first making /c/ontig /m/embe/s/hip /l/ists of each secondary cluster
#--'make.FASTA.file.set' is a function from our package RKXM (https://github.com/rbhwilliams/RKXM)
a.metadata.2_cmsl<-tapply(rownames(a.metadata.2),INDEX=a.metadata.2$sc.id,FUN=function(x){x})
make.FASTA.file.set(a.metadata.2_cmsl,a,"./a.metadata.2_cmsl_fna")

#--the example ends here--


