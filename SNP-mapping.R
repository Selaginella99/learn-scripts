setwd('~/Desktop/')

library(Hmisc)

## load gene lists to be mapped for their snps
sheet1 = read.csv('sheet1-aging.csv', as.is = TRUE)
sheet2 = read.csv('sheet2-longevity.csv', as.is = TRUE)
sheet3 = read.csv('sheet3-longevity.csv', as.is = TRUE)

# describe(mydata) to do some exploratory analysis

## load snp-gene mapping table
mapping = read.csv('IBC_chip_final_49K SNP_list_BroadInst.csv', as.is = TRUE)

## get uniqued gene symbol vector from sheet 1
gene1 = unique(sheet1$HGNC.Symbol)

## get uniqued gene symbol vector from sheet 2
## TODO: deal with the 'mutiple' ones
gene2 = unique(sheet2$Gene.Variant)[-which(unique(sheet2$Gene.Variant) == '')]

## get uniqued gene symbol vector from sheet 3
gene3 = unique(unlist(strsplit(sheet3$gene.s., ',')))

## create a list named genemap to store the known snp-gene mapping
## each snp in the data frame mapping may have several corresponding genes
## each component of the list stores a gene symbol vector for the snp
genemap = vector('list', nrow(mapping))
## split (and remove) the original gene symbol vector for comma
for (i in 1L:length(genemap)) genemap[[i]] = unlist(strsplit(mapping$RefSeq_hg18_Gene[i], ','))
## split (and remove) the original gene symbol vector for semicolon
for (i in 1L:length(genemap)) genemap[[i]] = unlist(strsplit(genemap[[i]], ';'))
## replace "NEAR GENE:" to "" for the original gene symbol vector
for (i in 1L:length(genemap)) genemap[[i]] = gsub('NEAR GENE:', '', genemap[[i]])
## replace "IN GENE:" to "" for the original gene symbol vector
for (i in 1L:length(genemap)) genemap[[i]] = gsub('IN GENE:', '', genemap[[i]])
## get the uniqued gene symbols for each snp
for (i in 1L:length(genemap)) genemap[[i]] = unique(genemap[[i]])

## get the uniqued gene symbol vector from sheet 1
candi_genes1 = unique(gene1)
## create a list named candi_snps1 to store the gene-snp mapping to be made
candi_snps1 = vector('list', length(candi_genes1))
## name each component of the list by the gene symbol
names(candi_snps1) = candi_genes1

## looking for the corresponding snps for each gene
for (i in 1L:length(candi_genes1)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes1[i] %in% genemap[[j]]) 
      candi_snps1[[i]] = c(candi_snps1[[i]], mapping$SNP_ID[j])
  }
}

## get the uniqued snp vector for sheet 1
candi_snps1 = unique(unlist(candi_snps1))

## get the uniqued gene symbol vector from sheet 2
candi_genes2 = unique(gene2)
## create a list named candi_snps2 to store the gene-snp mapping to be made
candi_snps2 = vector('list', length(candi_genes2))
## name each component of the list by the gene symbol
names(candi_snps2) = candi_genes2

## looking for the corresponding snps for each gene
for (i in 1L:length(candi_genes2)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes2[i] %in% genemap[[j]]) 
      candi_snps2[[i]] = c(candi_snps2[[i]], mapping$SNP_ID[j])
  }
}

## get the uniqued snp vector for sheet 2
candi_snps2 = unique(unlist(candi_snps2))

## get the uniqued gene symbol vector from sheet 3
candi_genes3 = unique(gene3)
## create a list named candi_snps3 to store the gene-snp mapping to be made
candi_snps3 = vector('list', length(candi_genes3))
## name each component of the list by the gene symbol
names(candi_snps3) = candi_genes3

## looking for the corresponding snps for each gene
for (i in 1L:length(candi_genes3)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes3[i] %in% genemap[[j]]) 
      candi_snps3[[i]] = c(candi_snps3[[i]], mapping$SNP_ID[j])
  }
}

## get the uniqued snp vector for sheet 3
candi_snps3 = unique(unlist(candi_snps3))

## get the uniqued gene symbol vector after combining sheet 1, 2, and 3
candi_genes = unique(c(gene1, gene2, gene3))
## create a list named candi_snps to store the gene-snp mapping to be made
candi_snps = vector('list', length(candi_genes))
## name each component of the list by the gene symbol
names(candi_snps) = candi_genes

## looking for the corresponding snps for each gene
for (i in 1L:length(candi_genes)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes[i] %in% genemap[[j]]) 
      candi_snps[[i]] = c(candi_snps[[i]], mapping$SNP_ID[j])
  }
}

## get the uniqued snp vector for all sheets
candi_snps = unique(unlist(candi_snps))

## write candidate snps, without gene names
write(candi_snps1, file = 'candi-snps-ibc-1.csv')
write(candi_snps2, file = 'candi-snps-ibc-2.csv')
write(candi_snps3, file = 'candi-snps-ibc-3.csv')
write(candi_snps, file = 'candi-snps-ibc-all.csv')

## write candidate snps, with gene names after
names(genemap) = mapping$SNP_ID
for (i in 1L:length(candi_snps1)) write(paste(candi_snps1[i], paste(genemap[[candi_snps1[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-1.csv', append = TRUE)
for (i in 1L:length(candi_snps2)) write(paste(candi_snps2[i], paste(genemap[[candi_snps2[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-2.csv', append = TRUE)
for (i in 1L:length(candi_snps3)) write(paste(candi_snps3[i], paste(genemap[[candi_snps3[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-3.csv', append = TRUE)
for (i in 1L:length(candi_snps))  write(paste(candi_snps[i], paste(genemap[[candi_snps[i]]], collapse = ','), sep = ','),   file = 'candi-snps-ibc-with-gene-symbol-all.csv', append = TRUE)
