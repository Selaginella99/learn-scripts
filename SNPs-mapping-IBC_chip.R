library('xlsx')
library(Hmisc)

setwd("./Candidate genes")

## load gene lists to be mapped for their snps

sheet1 = read.csv("GenAge.csv",as.is=TRUE)
sheet2 = read.csv("Longevity-link.csv",as.is=TRUE)
sheet3 = read.csv("LongevityMap.csv",as.is=TRUE)

## load snp-gene mapping table
mapping = read.csv('IBC_chip_final_49K SNP_list_BroadInst.csv', as.is = TRUE)

# describe(mydata) like describe(mapping) to explore the datasheet

candi_snp_1 = sheet1
for (i in 7:9) candi_snp_1[, i] = NA
names(candi_snp_1)[7:9] = c('SNP_ID', 'SNP_ID (in gene)', 'SNP_ID (near gene)')

candi_snp_2 = sheet2
for (i in 7:9) candi_snp_2[, i] = NA
names(candi_snp_2)[7:9] = c('SNP_ID', 'SNP_ID (in gene)', 'SNP_ID (near gene)')

candi_snp_3 = sheet3
for (i in 8:10) candi_snp_3[, i] = NA
names(candi_snp_3)[8:10] = c('SNP_ID', 'SNP_ID (in gene)', 'SNP_ID (near gene)')

## create three lists to store the known snp-gene mapping
## each snp in the data frame mapping may have several corresponding genes
## each component of the list stores a gene symbol vector for the snp
genemap_regular = vector('list', nrow(mapping))
genemap_ingene = vector('list', nrow(mapping))
genemap_neargene = vector('list', nrow(mapping))

for (i in 1:nrow(mapping)) {
  
  if (!grepl('NEAR GENE', mapping$RefSeq_hg18_Gene[i]) & 
        !grepl('IN GENE', mapping$RefSeq_hg18_Gene[i])) {
    genemap_regular[[i]] = mapping$RefSeq_hg18_Gene[i]
    genemap_ingene[[i]] = NA
    genemap_neargene[[i]] = NA
  }
  
  if (grepl('NEAR GENE', mapping$RefSeq_hg18_Gene[i]) &
        !grepl('IN GENE', mapping$RefSeq_hg18_Gene[i])) {
    genemap_regular[[i]] = NA
    genemap_ingene[[i]] = NA
    genemap_neargene[[i]] = unlist(strsplit(mapping$RefSeq_hg18_Gene[i], ','))
    genemap_neargene[[i]] = unlist(strsplit(genemap_neargene[[i]], ';'))
    genemap_neargene[[i]] = gsub('NEAR GENE:', '', genemap_neargene[[i]])
    genemap_neargene[[i]] = unique(genemap_neargene[[i]])
  }
  
  if (!grepl('NEAR GENE', mapping$RefSeq_hg18_Gene[i]) &
        grepl('IN GENE', mapping$RefSeq_hg18_Gene[i])) {
    genemap_regular[[i]] = NA
    genemap_neargene[[i]] = NA
    genemap_ingene[[i]] = unlist(strsplit(mapping$RefSeq_hg18_Gene[i], ','))
    genemap_ingene[[i]] = unlist(strsplit(genemap_ingene[[i]], ';'))
    genemap_ingene[[i]] = gsub('IN GENE:', '', genemap_ingene[[i]])
    genemap_ingene[[i]] = unique(genemap_ingene[[i]])
  }
  
  if (grepl('NEAR GENE', mapping$RefSeq_hg18_Gene[i]) &
        grepl('IN GENE', mapping$RefSeq_hg18_Gene[i])) {
    
    genemap_regular[[i]] = NA
    
    tmp0 = regexpr('(?<=IN GENE:)(.*?)+(?=NEAR GENE)', text = mapping$RefSeq_hg18_Gene[i], perl = TRUE)
    tmp1 = substring(mapping$RefSeq_hg18_Gene[i], tmp0, tmp0 + attr(tmp0, "match.length") - 1)
    genemap_ingene[[i]] = unique(unlist(strsplit(tmp1, ',')))
    genemap_ingene[[i]] = unlist(strsplit(genemap_ingene[[i]], ';'))
    genemap_ingene[[i]] = unique(genemap_ingene[[i]])
    
    tmp0 = regexpr('(?<=NEAR GENE:)(.*)+(?=;)', text = mapping$RefSeq_hg18_Gene[i], perl = TRUE)
    tmp1 = substring(mapping$RefSeq_hg18_Gene[i], tmp0, tmp0 + attr(tmp0, "match.length") - 1)
    genemap_neargene[[i]] = unique(unlist(strsplit(tmp1, ',')))
    genemap_neargene[[i]] = unlist(strsplit(genemap_neargene[[i]], ';'))
    genemap_neargene[[i]] = unique(genemap_neargene[[i]])
    
  }
  
}

## looking for the corresponding snps for each gene
for (i in 1L:nrow(candi_snp_1)) {
  print(i)
  for (j in 1L:length(genemap_regular)) {
    if (candi_snp_1$HGNC.Symbol[i] %in% genemap_regular[[j]])
      candi_snp_1[i, 'SNP_ID'] = paste(na.omit(c(candi_snp_1[i, 'SNP_ID'], mapping$SNP_ID[j])), collapse = ',')
    
    if (candi_snp_1$HGNC.Symbol[i] %in% genemap_ingene[[j]])
      candi_snp_1[i, 'SNP_ID (in gene)'] = paste(na.omit(c(candi_snp_1[i, 'SNP_ID (in gene)'], mapping$SNP_ID[j])), collapse = ',')
    
    if (candi_snp_1$HGNC.Symbol[i] %in% genemap_neargene[[j]])
      candi_snp_1[i, 'SNP_ID (near gene)'] = paste(na.omit(c(candi_snp_1[i, 'SNP_ID (near gene)'], mapping$SNP_ID[j])), collapse = ',')
  }
}

## looking for the corresponding snps for each gene
for (i in 1L:nrow(candi_snp_2)) {
  print(i)
  for (j in 1L:length(genemap_regular)) {
    if (candi_snp_2$Gene.Variant[i] %in% genemap_regular[[j]])
      candi_snp_2[i, 'SNP_ID'] = paste(na.omit(c(candi_snp_2[i, 'SNP_ID'], mapping$SNP_ID[j])), collapse = ',')
    
    if (candi_snp_2$Gene.Variant[i] %in% genemap_ingene[[j]])
      candi_snp_2[i, 'SNP_ID (in gene)'] = paste(na.omit(c(candi_snp_2[i, 'SNP_ID (in gene)'], mapping$SNP_ID[j])), collapse = ',')
    
    if (candi_snp_2$Gene.Variant[i] %in% genemap_neargene[[j]])
      candi_snp_2[i, 'SNP_ID (near gene)'] = paste(na.omit(c(candi_snp_2[i, 'SNP_ID (near gene)'], mapping$SNP_ID[j])), collapse = ',')
  }
}

## looking for the corresponding snps for each gene
for (i in 1L:nrow(candi_snp_3)) {
  
  print(i)
  
  gene_names = unlist(strsplit(candi_snp_3$gene.s.[i], ','))
  
  if (length(gene_names) > 1L) {
    
    for (j in 1L:length(genemap_regular)) {
      for (k in 1L:length(gene_names)) {
        if (gene_names[k] %in% genemap_regular[[j]])
          candi_snp_3[i, 'SNP_ID'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID'], mapping$SNP_ID[j])), collapse = ',')
        
        if (gene_names[k] %in% genemap_ingene[[j]])
          candi_snp_3[i, 'SNP_ID (in gene)'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID (in gene)'], mapping$SNP_ID[j])), collapse = ',')
        
        if (gene_names[k] %in% genemap_neargene[[j]]) 
          candi_snp_3[i, 'SNP_ID (near gene)'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID (near gene)'], mapping$SNP_ID[j])), collapse = ',')
      }
    }
    
    if (!is.na(candi_snp_3[i, 'SNP_ID'])) candi_snp_3[i, 'SNP_ID'] = paste(unique(unlist(strsplit(candi_snp_3[i, 'SNP_ID'], ','))), collapse = ',')
    if (!is.na(candi_snp_3[i, 'SNP_ID (in gene)'])) candi_snp_3[i, 'SNP_ID (in gene)'] = paste(unique(unlist(strsplit(candi_snp_3[i, 'SNP_ID (in gene)'], ','))), collapse = ',')
    if (!is.na(candi_snp_3[i, 'SNP_ID (near gene)'])) candi_snp_3[i, 'SNP_ID (near gene)'] = paste(unique(unlist(strsplit(candi_snp_3[i, 'SNP_ID (near gene)'], ','))), collapse = ',')
    
  } else {
    
    for (j in 1L:length(genemap_regular)) {
      if (candi_snp_3$gene.s.[i] %in% genemap_regular[[j]])
        candi_snp_3[i, 'SNP_ID'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID'], mapping$SNP_ID[j])), collapse = ',')
      
      if (candi_snp_3$gene.s.[i] %in% genemap_ingene[[j]])
        candi_snp_3[i, 'SNP_ID (in gene)'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID (in gene)'], mapping$SNP_ID[j])), collapse = ',')
      
      if (candi_snp_3$gene.s.[i] %in% genemap_neargene[[j]])
        candi_snp_3[i, 'SNP_ID (near gene)'] = paste(na.omit(c(candi_snp_3[i, 'SNP_ID (near gene)'], mapping$SNP_ID[j])), collapse = ',')
    }
    
  }
  
}

## write candidate snps in xlsx format
write.xlsx(candi_snp_1, file = 'candi-snps-ibc-1.xlsx', showNA = FALSE, row.names = FALSE)
write.xlsx(candi_snp_2, file = 'candi-snps-ibc-2.xlsx', showNA = FALSE, row.names = FALSE)
write.xlsx(candi_snp_3, file = 'candi-snps-ibc-3.xlsx', showNA = FALSE, row.names = FALSE)

## names(genemap) = mapping$SNP_ID
# for (i in 1L:length(candi_snps1)) write(paste(candi_snps1[i], paste(genemap[[candi_snps1[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-1.csv', append = TRUE)
# for (i in 1L:length(candi_snps2)) write(paste(candi_snps2[i], paste(genemap[[candi_snps2[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-2.csv', append = TRUE)
# for (i in 1L:length(candi_snps3)) write(paste(candi_snps3[i], paste(genemap[[candi_snps3[i]]], collapse = ','), sep = ','), file = 'candi-snps-ibc-with-gene-symbol-3.csv', append = TRUE)
# for (i in 1L:length(candi_snps))  write(paste(candi_snps[i], paste(genemap[[candi_snps[i]]], collapse = ','), sep = ','),   file = 'candi-snps-ibc-with-gene-symbol-all.csv', append = TRUE)


