genoprefix = 'phg000004.ind.geno.'  # set the prefix of the genotype files here

famfile = list.files()[grep('.fam', list.files())]
famtable = read.table(famfile, sep = ' ', header = FALSE, as.is = TRUE)
indindex = as.character(famtable$V2)

snps = read.table(paste0(paste0(genoprefix, indindex[1L]), '.gz'), sep = ',', header = FALSE, as.is = TRUE)[, 1L]

evoker_df = as.data.frame(matrix(NA, nrow = length(snps), ncol = (2L * length(indindex)) + 1L))
names(evoker_df) = c('SNP', rep(indindex, each = 2L))
evoker_df$SNP = snps

for (i in 1L:length(indindex)) {
  cat('Merging', i, 'in', length(indindex), '\n')
  tmp = read.table(paste0(paste0(genoprefix, indindex[i]), '.gz'),
                   sep = ',', header = FALSE, 
                   colClasses = c(rep('NULL', 8L), 
                                  rep('numeric', 2)))  # load faster
  evoker_df[, 2L * i] = tmp[, 1L]
  evoker_df[, 2L * i + 1L] = tmp[, 2L]
}

write.table(evoker_df, file = 'input.txt', sep = ' ', quote = FALSE, col.names = TRUE, row.names = FALSE)
