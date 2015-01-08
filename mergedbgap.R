setwd('~/Downloads/Test_Data/')
filenames = list.files()
objnames = gsub('_', '.', gsub('-', '.', gsub('.DBGAP.TXT.ind', '', filenames)))

for (i in 1L:length(filenames)) {
  assign(objnames[i], read.table(filenames[i], sep = ' ', as.is = TRUE, 
                                 colClasses = c(rep('character', 3L), rep('NULL', 3L))))
}

x = eval(parse(text = paste0("cbind(", paste(paste0(objnames, '[, 2:3]'), collapse = ', '), ")")))
x = eval(parse(text = paste0("cbind(", objnames[1L], "[, 1], x)")))
names(x) = c('SNP.id', paste(rep(objnames, each = 2L), c('allele1', 'allele2'), sep = '.'))
