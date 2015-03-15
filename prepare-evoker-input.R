genoprefix = 'phg000004.ind.geno.'  # set the prefix of the genotype files here

famfile = list.files()[grep('.fam', list.files())]
famtable = read.table(famfile[1], sep = ' ', header = FALSE, as.is = TRUE)
indindex = as.character(famtable$V2)
# indindex = c('5', '9', '26', '33', '35')  # for testing

snps = read.table(paste0(paste0(genoprefix, indindex[1L]), '.gz'), sep = ',', header = FALSE, as.is = TRUE)[, 1L]

evoker_df = as.data.frame(matrix(NA, nrow = length(snps), ncol = (2L * length(indindex)) + 1L))
names(evoker_df) = c('SNP', rep(indindex, each = 2L))
evoker_df$SNP = snps

for (i in 1L:length(indindex)) {
   cat('Loading', i, 'in', length(indindex), '\n')
   eval(parse(text = paste0('tmp_', i, ' = read.table(paste0(paste0(genoprefix, indindex[', i, ']), ".gz"), 
   sep = ",", header = FALSE, colClasses = c(rep("NULL", 8L), rep("character", 2L)))')))
 }
 
gc()
evoker_mat = do.call(cbind, lapply(paste0('tmp_', 1L:length(indindex)), get))
evoker_mat = cbind(snps, evoker_mat)
colnames(evoker_mat) = c('SNP', rep(indindex, each = 2L))

write.table(evoker_mat, file = 'input.txt', sep = ' ', quote = FALSE, col.names = TRUE, row.names = FALSE)
