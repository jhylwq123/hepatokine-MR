library(readxl)
library(dplyr)
library(coloc)
setwd("/public/home/jiahanying/project/rare_common/coloc/gwas_pQTL")
#decode_symbol_mr <- read.delim("~/project/rare_common/coloc/gwas_pQTL/decode_symbol_mr.txt")
combine <- read.delim("~/project/rare_common/coloc/gwas_pQTL/GWAS_1393_r0.1_10000kb_100kb_pQTL_0.01_decode.txt")
combine[combine$P_pQTL==0,"Pval_pQTL"] <- 3.50e-320
combine[combine$P_GWAS==0,"P_GWAS"] <- 3.509842e-320
#combine <- merge(combine,decode_symbol_mr,by.x="Gene",by.y="Gene")
#write.table(combine,"GWAS_1393_r0.1_10000kb_250kb_pQTL_0.01_decode_mr.txt",sep="\t",quote=F,row.names=F)
list <- unique(combine$lead_snp_gene)
coloc=data.frame()
for (i in 1:length(list)){
    SNP_gene=list[i]
    data=combine[combine$lead_snp_gene==SNP_gene,]
    data <- data[!duplicated(data$SNP),]
    result=coloc.abf(dataset1 = list(pvalues=data$P_GWAS,snp=data$SNP,type="quant",N=1466352),
                     dataset2 = list(pvalues=data$Pval_pQTL,snp=data$SNP,type="quant",N=45000),MAF = data$MAF_GWAS)
    PPH4abf <- as.data.frame(result$summary)[6,1]
    PPH3abf <- as.data.frame(result$summary)[5,1]
    lead_snp_gene <- SNP_gene
    table <- result$results
    table$lead_snp_gene <- SNP_gene
    table$PPH4abf <- PPH4abf
    table$PPH3abf <- PPH3abf
    table$nsnp <-nrow(data)
    print(i)
    coloc<-rbind(coloc,table)
}
write.table(coloc,"coloc_gwas_0.1_10000kb_100kb_decode_0.01_all_full_results.txt",sep="\t",quote=F,row.names=F)
