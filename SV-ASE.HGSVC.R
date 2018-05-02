

# Begin SV-ASE analysis for integrated illumina IL-SV and merged PB-SV callsets.
# First read in the results from bedtools multicov.

#counts<-read.table(file="HG00514.illumina.hetDELINS.multicov.040318.txt", stringsAsFactors = FALSE)
#counts<-read.table(file="HG00733.illumina.hetDELINS.multicov.040318.txt", stringsAsFactors = FALSE)
#counts<-read.table(file="NA19240.illumina.hetDELINS.multicov.040318.txt", stringsAsFactors = FALSE)
#counts<-read.table(file="HG00514.INV.multicov..txt", stringsAsFactors = FALSE)
#counts<-read.table(file="HG00733.INV.multicov.txt", stringsAsFactors = FALSE)
counts<-read.table(file="NA19240.INV.multicov.txt", stringsAsFactors = FALSE)

colnames(counts)<-c("gene.chrom","gene.start","gene.stop","geneID","geneName","sv.chrom","sv.start","sv.end","ref","alt","qual","filt","sampleID","GT","hap1","hap2") 
counts$ratio<-apply(counts,1,function(x) min(as.numeric(x[15]),as.numeric(x[16]))/max(as.numeric(x[15]),as.numeric(x[16])))
counts$total<-apply(counts,1,function(x) sum(as.numeric(x[15]),as.numeric(x[16])))
counts$pval<-apply(counts,1,function(x) binom.test(as.numeric(x[15]),as.numeric(x[18]))$p.value)
myp<-counts$pval
counts$FDR<-p.adjust(myp,method="fdr")
SVsig<-counts[counts$FDR<=0.05,]
Sv<-SVsig[with(SVsig,order(as.numeric(FDR))),]
write.table(Sv,file="NA19240.illumina.hetINV.SVASE.txt",sep="\t",quote=F,col.names=TRUE,row.names=FALSE)

# Minor variation for format of the PB-SV callset shown below, the main analysis is identical.

#counts<-read.table(file="HG00514.PBSV.multicov.040318.txt", stringsAsFactors = FALSE)
#counts<-read.table(file="HG00733.PBSV.multicov.040318.txt", stringsAsFactors = FALSE)
counts<-read.table(file="NA19240.PBSV.multicov.040318.txt", stringsAsFactors = FALSE)

colnames(counts)<-c("gene.chrom","gene.start","gene.stop","geneID","geneName","sv.chrom","sv.start","sv.end","ref","alt","qual","filt","svLength","gt","GT","hap1","hap2") 
counts$ratio<-apply(counts,1,function(x) min(as.numeric(x[16]),as.numeric(x[17]))/max(as.numeric(x[16]),as.numeric(x[17])))
counts$total<-apply(counts,1,function(x) sum(as.numeric(x[16]),as.numeric(x[17])))
counts$pval<-apply(counts,1,function(x) binom.test(as.numeric(x[16]),as.numeric(x[19]))$p.value)
myp<-counts$pval
counts$FDR<-p.adjust(myp,method="fdr")
SVsig<-counts[counts$FDR<=0.05,]
Sv<-SVsig[with(SVsig,order(as.numeric(FDR))),]
write.table(Sv,file="NA19240.pacbio.SVASE.txt",sep="\t",quote=F,col.names=TRUE,row.names=FALSE)





