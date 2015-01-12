setwd('~/Desktop/')

sheet1 = read.csv('sheet1-aging.csv', as.is = TRUE)
sheet2 = read.csv('sheet2-longevity.csv', as.is = TRUE)
sheet3 = read.csv('sheet3-longevity.csv', as.is = TRUE)

mapping = read.csv('IBC_chip_final_49K SNP_list_BroadInst.csv', as.is = TRUE)

gene1 = unique(sheet1$HGNC.Symbol)

# TODO: deal with the 'mutiple' ones
gene2 = unique(sheet2$Gene.Variant)[-which(unique(sheet2$Gene.Variant) == '')]

gene3 = unique(unlist(strsplit(sheet3$gene.s., ',')))

##
genemap = vector('list', nrow(mapping))
for (i in 1L:length(genemap)) genemap[[i]] = unlist(strsplit(mapping$RefSeq_hg18_Gene[i], ','))
for (i in 1L:length(genemap)) genemap[[i]] = unlist(strsplit(genemap[[i]], ';'))
for (i in 1L:length(genemap)) genemap[[i]] = unlist(strsplit(genemap[[i]], ';'))
for (i in 1L:length(genemap)) genemap[[i]] = gsub('NEAR GENE:', '', genemap[[i]])
for (i in 1L:length(genemap)) genemap[[i]] = gsub('IN GENE:', '', genemap[[i]])
for (i in 1L:length(genemap)) genemap[[i]] = unique(genemap[[i]])

##
candi_genes1 = unique(gene1)
candi_snps1 = vector('list', length(candi_genes1))
names(candi_snps1) = candi_genes1

for (i in 1L:length(candi_genes1)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes1[i] %in% genemap[[j]]) 
      candi_snps1[[i]] = c(candi_snps1[[i]], mapping$SNP_ID[j])
  }
}

candi_snps1 = unique(unlist(candi_snps1))

##
candi_genes2 = unique(gene2)
candi_snps2 = vector('list', length(candi_genes2))
names(candi_snps2) = candi_genes2

for (i in 1L:length(candi_genes2)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes2[i] %in% genemap[[j]]) 
      candi_snps2[[i]] = c(candi_snps2[[i]], mapping$SNP_ID[j])
  }
}

candi_snps2 = unique(unlist(candi_snps2))

##
candi_genes3 = unique(gene3)
candi_snps3 = vector('list', length(candi_genes3))
names(candi_snps3) = candi_genes3

for (i in 1L:length(candi_genes3)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes3[i] %in% genemap[[j]]) 
      candi_snps3[[i]] = c(candi_snps3[[i]], mapping$SNP_ID[j])
  }
}

candi_snps3 = unique(unlist(candi_snps3))

##
candi_genes = unique(c(gene1, gene2, gene3))
candi_snps = vector('list', length(candi_genes))
names(candi_snps) = candi_genes

for (i in 1L:length(candi_genes)) {
  print(i)
  for (j in 1L:length(genemap)) {
    if (candi_genes[i] %in% genemap[[j]]) 
      candi_snps[[i]] = c(candi_snps[[i]], mapping$SNP_ID[j])
  }
}

candi_snps = unique(unlist(candi_snps))

write(candi_snps1, file = 'candi-snps-1.csv')
write(candi_snps2, file = 'candi-snps-2.csv')
write(candi_snps3, file = 'candi-snps-3.csv')
write(candi_snps, file = 'candi-snps.csv')
