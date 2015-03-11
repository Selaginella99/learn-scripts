## Find the overlap of SNPs between male & female for LS & RS and map back to genes
## Date: March 11, 2015

snp_ls_f = read.csv('Markers LS F_chi.csv', header = TRUE, as.is = TRUE)  # load snp of lifespan table of female (continuous outcome)
snp_ls_m = read.csv('Markers LS M_chi.csv', header = TRUE, as.is = TRUE)  # load snp of lifespan table of male (continuous outcome)
snp_rs_f = read.csv('Markers RS F_chi.csv', header = TRUE, as.is = TRUE)  # load snp table of female (continuous outcome)
snp_rs_m = read.csv('Markers RS M_chi.csv', header = TRUE, as.is = TRUE)  # load snp table of male (continuous outcome)

ibctable = read.csv('IBC_chip_final_49K SNP_list_BroadInst.csv', header = TRUE, as.is = TRUE)
row.names(ibctable) = ibctable$SNP_ID

marker_ls_f = snp_ls_f$Marker  # extract marker from table: ls female snp id (character vector)
marker_ls_m = snp_ls_m$Marker  # extract marker from table: ls male snp id (character vector)
marker_rs_f = snp_rs_f$Marker  # extract marker from table: rs female snp id (character vector)
marker_rs_m = snp_rs_m$Marker  # extract marker from table: rs male snp id (character vector)

ls_intersect = intersect(marker_ls_f, marker_ls_m)  # get intersect set of M & F for LS
rs_intersect = intersect(marker_rs_f, marker_rs_m)  # get intersect set of M & F for RS

ls_final_table = ibctable[ls_intersect, c('SNP_ID', 'Chromosome', 'RefSeq_hg18_Gene')]
rs_final_table = ibctable[rs_intersect, c('SNP_ID', 'Chromosome', 'RefSeq_hg18_Gene')]

ls_na_id = which(is.na(ls_final_table$SNP_ID))  
rs_na_id = which(is.na(rs_final_table$SNP_ID))  
ls_final_table[ls_na_id, 'SNP_ID'] = ls_intersect[ls_na_id]  # refill snp id to the NA rows of the snp id column
rs_final_table[rs_na_id, 'SNP_ID'] = rs_intersect[rs_na_id]  # refill snp id to the NA rows of the snp id column

write.csv(ls_final_table, file = 'intersect_table_ls.csv', row.names = FALSE, quote = FALSE)
write.csv(rs_final_table, file = 'intersect_table_rs.csv', row.names = FALSE, quote = FALSE)
