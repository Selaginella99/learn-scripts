## Mapping the intersecting SNPs/genes in the EMMAX to the annotated IBC chip
library(Hmisc)

setwd("Z:/Project/Framingham/Data/Phenotype/FHS49K_lifespan&weight")

## load snps/gene lists to be mapped for their SNPs/genes in the annotation of IBC49K chip

lr_e_f = read.csv("LR_emmax_F.csv",as.is=TRUE)
lr_e_m = read.csv("LR_emmax_M.csv",as.is=TRUE)

## load snp-gene mapping table
mapping.ibc = read.csv('IBC_chip_final_49K SNP_list_BroadInst.csv', as.is = TRUE)
# describe(mydata) like describe(mapping) to explore the datasheet

## set up the new snp-gene mapping table
intersect_emmax_gene_1 = lr_e_f
for (i in 5) intersect_emmax_gene_1[, i] = NA
names(intersect_emmax_gene_1)[5] = c('IBC_gene')

### intersect_emmax_gene_1 = lr_e_f
### intersect_emmax_gene_1[, 5] = NA
### names(intersect_emmax_gene_1)[5] = c('IBC_gene')

intersect_emmax_gene_2 = lr_e_m
for (i in 5) intersect_emmax_gene_2[, i] = NA
names(intersect_emmax_gene_2)[5] = c('IBC_gene')

### intersect_emmax_gene_2 = lr_e_m
### intersect_emmax_gene_2[, 5] = NA
### names(intersect_emmax_gene_2)[5] = c('IBC_gene')


## looking for the corresponding gene for each SNP
for (i in 1L:nrow(intersect_emmax_gene_1)) {
  print(i)
  for (j in 1L:nrow(mapping.ibc)) {
    if (intersect_emmax_gene_1$rs_num[i] %in% mapping.ibc$SNP_ID[j])
      intersect_emmax_gene_1[i, 'IBC_gene'] = paste(na.omit(c(intersect_emmax_gene_1[i, 'IBC_gene'], mapping.ibc$RefSeq_hg18_Gene[j])), collapse = ',')    
  }
}

for (i in 1L:nrow(intersect_emmax_gene_2)) {
  print(i)
  for (j in 1L:nrow(mapping.ibc)) {
    if (intersect_emmax_gene_2$rs_num[i] %in% mapping.ibc$SNP_ID[j])
      intersect_emmax_gene_2[i, 'IBC_gene'] = paste(na.omit(c(intersect_emmax_gene_2[i, 'IBC_gene'], mapping.ibc$RefSeq_hg18_Gene[j])), collapse = ',')    
  }
}

## write candidate snps in csv format
write.csv(intersect_emmax_gene_1, file = 'LR_Emmax_F_ibc.csv')
write.csv(intersect_emmax_gene_2, file = 'LR_Emmax_M_ibc.csv')
