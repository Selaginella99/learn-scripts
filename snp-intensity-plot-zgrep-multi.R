setwd('~/Desktop/fhs/')  # directory that stores the genotype files
ssid = c('ss74822184', 'ss74805882', 'ss74820711', 'ss74814609')  # set ssid of snps to plot
rsid = c('rs2233595', 'rs9928967', 'rs8081943', 'rs1794108')  # set corresponding rsid (must be coherent with ssid)
genofilepattern = 'phg000004.ind.geno.*?.gz'  # set the file name pattern of genotype files

# get the names of gzipped genotype files
genofiles = list.files()[grep(genofilepattern, list.files())]

# initialize the list for storing the intensity data
intlist = vector('list', length(rsid))
names(intlist) = rsid
for (i in 1:length(intlist)) intlist[[i]] = as.data.frame(matrix(NA, ncol = 4, nrow = length(genofiles)))
for (i in 1:length(intlist)) names(intlist[[i]]) = c('x', 'y', 'Genotype', 'individual')

# extract the intensity data from genotype files
# non-parallelized, but zgrep multiple files + multiple snps at the same time
pattern = paste(paste(ssid, rsid, sep = ','), collapse = '|')
doc = system(paste0("LC_ALL=C find . -name '", genofilepattern, 
                    "' -exec zgrep -H -E '", pattern, "' {} \\;"), intern = TRUE)
# example: LC_ALL=C find . -name 'phg000004.ind.geno.*?.gz' -exec zgrep -H -E 'ss74822184,rs2233595|ss74805882,rs9928967|ss74820711,rs8081943' {} \;

# tidy result
txt = gsub('./phg000004.ind.geno.|.gz', '', gsub(':', ',', doc))

# load result into list
for (k in 1:length(rsid)) {
  tmp = as.matrix(data.frame(strsplit(txt[grep(paste(ssid[k], rsid[k], sep = ','), txt)], ',')))
  colnames(tmp) = NULL
  intlist[[rsid[k]]][, 1] = as.numeric(tmp[10, ])  # x
  intlist[[rsid[k]]][, 2] = as.numeric(tmp[11, ])  # y
  intlist[[rsid[k]]][, 3] = tmp[6, ]  # genotype
  intlist[[rsid[k]]][, 4] = tmp[1, ]  # a clean individual id
}

# save result to hdd
save(intlist, file = 'intlist.RData')

library('ggplot2')
library('ellipse')

# function for plotting the intensity data
plot.intensity = function(rsid) {
  
  intdf = intlist[[rsid]]
  
  if (length(which(intdf$Genotype == 'NN')) != 0L) {
    
    # set 'NN' to '--'
    intdf[which(intdf$Genotype == 'NN'), 'Genotype'] = '--'
    
    # move NN to the last position
    genofactor = factor(intdf$Genotype)
    nnlocation = which(levels(genofactor) == '--')
    intdf$Genotype = factor(intdf$Genotype, 
                            levels = levels(genofactor)[c((1:nlevels(genofactor))[-nnlocation], nnlocation)], 
                            ordered = TRUE)
    
  }
  
  # define custom legend labels
  genolabel = paste0(levels(intdf$Genotype), ' (', summary(intdf$Genotype), ')')
  
  # generating ellipses
  df_ell = data.frame()
  for(g in levels(as.factor(intdf$Genotype))) {
    df_ell = rbind(df_ell, 
                   cbind(as.data.frame(with(intdf[intdf$Genotype == g, ], 
                                            ellipse(cor(x, y), 
                                                    scale = c(sd(x), sd(y)), 
                                                    centre = c(mean(x), mean(y))))), 
                         Genotype = g))
  }
  
  pal = c('#d62728','#e377c2','#1f77b4','#7f7f7f')
  # pal = c('#1f77b4', '#ff7f0e', '#2ca02c', '#7f7f7f')
  # color palette from: https://github.com/mbostock/d3/wiki/Ordinal-Scales
  
  p = ggplot(data = intdf, aes(x = x, y = y, colour = Genotype)) +
    geom_point(size = 1, alpha = 0.7, show_guide = FALSE) +
    geom_path(data = df_ell, aes(x = x, y = y, colour = Genotype), size = 0.5, linetype = 2) +
    scale_colour_manual(values = pal, labels = genolabel) +
    labs(title = paste0(rsid, ', All individuals'),
         x = 'Allele A (C) Intensity',
         y = 'Allele B (G) Intensity') +
    theme_bw() +
    theme(legend.key = element_blank())
  
  print(p)
  ggsave(paste0(rsid, '.pdf'), p, width = 10, height = 8)
  
}

# examples
plot.intensity(rsid[1])
plot.intensity(rsid[2])
plot.intensity(rsid[3])
plot.intensity(rsid[4])
