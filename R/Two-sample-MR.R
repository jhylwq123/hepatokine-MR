library(readxl)
library(TwoSampleMR)
library(MRPRESSO)
ukbb_1870_pQTL_GWAS_no_drug_target_gene_0.1_10000kb_summary <- read.delim("all_gene_gwas_pQTL_r2_0.1_10000kb.txt")
gene_samplesize <- read.delim("~/project/rare_common/two-sample-mr/blood_pQTL/input/clump_by_gene/gene_samplesize.txt")
pQTL_GWAS_no_drug_target_gene <- merge(ukbb_1870_pQTL_GWAS_no_drug_target_gene_0.1_10000kb_summary,gene_samplesize,by.x = "GENE",by.y = "GENE")
##formatting
#exp
exp_r2_0.1_pQTL <- pQTL_GWAS_no_drug_target_gene[,c(3,16,17,19,20,25,18,24,26)]
colnames(exp_r2_0.1_pQTL) <- c("SNP","other_allele.exposure","effect_allele.exposure","beta.exposure",
                                        "se.exposure","pval.exposure","eaf.exposure","gene_symbol","samplesize.exposure")
exp_r2_0.1_pQTL$exposure <- "pQTL"
exp_r2_0.1_pQTL$id.exposure <- "1"
#out
out_r2_0.1_pQTL <- pQTL_GWAS_no_drug_target_gene[c(3,6,7,9,11,12,13,14,2,15)]
colnames(out_r2_0.1_pQTL) <- c("SNP","other_allele.outcome","effect_allele.outcome","eaf.outcome","beta.outcome",
                                       "se.outcome","pval.outcome","samplesize.outcome","hg37","hg38")
out_r2_0.1_pQTL$outcome <- "eGFR"
out_r2_0.1_pQTL$id.outcome <- "2"
exp <- exp_r2_0.1_pQTL
out <- out_r2_0.1_pQTL
out <- out[!duplicated(out$SNP),]
##gene position
gene_coordinate_hg38 <- read.delim("~/project/rare_common/two-sample-mr/cis_blood_pQTL/input/gene_coordinate_hg38.txt")
##
drug_list <- exp[!duplicated(exp$gene_symbol),]
drug_list <- as.data.frame(drug_list)
##if
pvalu <-c()
me <- c()
drug <- c()
nsnp <- c()
beta <- c()
se <- c()
Egger_intercept <- c()
Egger_Q <- c()
Cochran_Q <- c()
for (i in 1:nrow(drug_list)){
  drug_id <- drug_list[i,8]
  exp_drug <- exp[exp$gene_symbol==drug_id,]
  exp_drug <- exp_drug[!duplicated(exp_drug$SNP),] ##avoid duplicated gene
  dat <- harmonise_data(exposure_dat = exp_drug, outcome_dat = out)
  dat <- merge(dat,gene_coordinate_hg38,by.x = "gene_symbol",by.y = "symbol") 
  dat[c("chrom", "bp")] <- as.data.frame(do.call(rbind, strsplit(as.character(dat[, 19]), "_", fixed = TRUE))) 
  dat[,34] <-as.numeric(as.character(dat[, 34]))
  dat_cis <- dat[dat[, 34] >= dat[, 31]-1000000 & dat[, 34] <= dat[, 32]+1000000, ] ##choose cis-QTL
  test <- dat_cis[dat_cis$mr_keep=="TRUE",]
  if (nrow(test)>= 4){
  res <- mr(dat_cis) ##MR分析
  pleio <- mr_pleiotropy_test(dat_cis) 
  het <- mr_heterogeneity(dat_cis)
  print(i)
  res$drug <- drug_list[i,8]
  res$Egger_intercept<- pleio[1,7]
  res$Egger_Q <- het[1,8]
  res$Cochran_Q <- het[2,8]
  p.val=res[,9]
  method=res[,5]
  dr=res[,10]
  n=res[,6]
  beta.1 =res[,7]
  se.1 =res[,8]
  Egger_intercept.1 = res[,11]
  Egger_Q.1 = res[,12]
  Cochran_Q.1 = res[,13]
   drug = c(drug,dr)
  me = c(me,method)
  pvalu = c(pvalu,p.val)
  beta = c(beta,beta.1)
  se = c(se,se.1)
  nsnp = c(nsnp,n)
  Egger_intercept = c(Egger_intercept,Egger_intercept.1)
  Egger_Q =c(Egger_Q,Egger_Q.1)
  Cochran_Q =c(Cochran_Q,Cochran_Q.1)
}
}
out_pQTL_drug <- data.frame(pvalu=pvalu,me=me,drug=drug,nsnp=nsnp,beta=beta,se=se,Egger_intercept=Egger_intercept,Egger_Q=Egger_Q,Cochran_Q=Cochran_Q)
write.table(out_pQTL_drug,"ukbb_cis_out_pQTL_r2_0.1_gwas_pQTL_above3.txt",sep="\t",quote=F,row.names=F)
