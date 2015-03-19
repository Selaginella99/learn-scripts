library('ellipse')
library('ggplot2')

rsid = 'rs00000000'

set.seed(101)
x = rnorm(500, mean = 2)
y = 1.5 + 0.4 * x + rnorm(500)
df = data.frame(x = x, y = y, Genotype = paste0("CC (", length(x), ')'))
x = rnorm(300, mean = 2)
y = 1.5 * x + 0.4 + rnorm(300)
df = rbind(df, data.frame(x = x, y = y, Genotype = paste0("CG (", length(x), ')')))
x = rnorm(200, mean = 2)
y = 4 * x + 2 + rnorm(200)
df = rbind(df, data.frame(x = x, y = y, Genotype = paste0("GG (", length(x), ')')))
x = rnorm(50, mean = 4)
y = rnorm(50)
df = rbind(df, data.frame(x = x, y = y, Genotype = paste0("-- (", length(x), ')')))

pal = c('#1f77b4', '#ff7f0e', '#2ca02c', '#7f7f7f')
# From: https://github.com/mbostock/d3/wiki/Ordinal-Scales

# Generating ellipses
df_ell = data.frame()
for(g in levels(df$Genotype)) {
  df_ell = rbind(df_ell, 
                 cbind(as.data.frame(with(df[df$Genotype == g, ], 
                                          ellipse(cor(x, y), 
                                                  scale = c(sd(x),sd(y)), 
                                                  centre = c(mean(x), mean(y))))), 
                       Genotype = g))
}

p = ggplot(data = df, aes(x = x, y = y, colour = Genotype)) +
  geom_point(size = 1.5, alpha = 0.7, show_guide = FALSE) +
  geom_path(data = df_ell, aes(x = x, y = y, colour = Genotype), size = 1, linetype = 2) +
  scale_colour_manual(values = pal) +
  labs(title = paste0(rsid, ', Original Cohort'),
       x = 'Allele A (C) Intensity',
       y = 'Allele B (G) Intensity') +
  theme_bw() +
  theme(legend.key = element_blank())

print(p)
ggsave(paste0('intensity-', rsid, '.pdf'), p, width = 10, height = 8)

##################################################################################################################

### read data and plot the SNP cluster plots
library('ggplot2')
library('ellipse')

setwd('~/Desktop/fhs/')  # directory that stores the genotype files
ssid = c("ss74814609","ss74820711","ss74805882","ss74822184") # set ssid of snps to plot
rsid = c('rs1794108', 'rs8081943', 'rs9928967','rs2233595')  # set corresponding rsid (must be coherent with ssid)
genoprefix = 'phg000004.ind.geno.'  # set the prefix of genotype files

genofiles = list.files()[grep(genoprefix, list.files())]  # please only keep .gz files in the directory

# initialize the list for storing the intensity data
intlist = vector('list', length(rsid))
names(intlist) = rsid
for (i in 1:length(intlist)) intlist[[i]] = as.data.frame(matrix(NA, ncol = 3, nrow = length(genofiles)))
for (i in 1:length(intlist)) names(intlist[[i]]) = c('x', 'y', 'Genotype')

# extract the intensity data from genotype files
for (j in 1:length(genofiles)) {
  cat('Reading', j, 'in', length(genofiles), '\n')
  tmp = read.table(genofiles[j], sep = ',', header = FALSE,
                   colClasses = c('character', 'character', 'NULL', 'NULL',
                                  'character', 'NULL', 'NULL', 'NULL',
                                  'numeric', 'numeric'))
  snpuniqueid = paste(tmp[, 1], tmp[, 2], sep = '_')
  row.names(tmp) = snpuniqueid
  for (k in 1:length(rsid)) {
    intlist[[rsid[k]]][j, ] = tmp[paste(ssid[k], rsid[k], sep = '_'), c(4, 5, 3)]
  }
}

# save the result to hdd
save(intlist, file = 'intlist.RData')

# function for plotting the intensity data
plot.intensity = function(rsid) {
  
  intdf = intlist[[rsid]]
  
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
  #pal = c('#1f77b4', '#ff7f0e', '#2ca02c', '#7f7f7f')
  # color palette from: https://github.com/mbostock/d3/wiki/Ordinal-Scales
  
  p = ggplot(data = intdf, aes(x = x, y = y, colour = Genotype)) +
    geom_point(size = 1.5, alpha = 0.7, show_guide = FALSE) +
    geom_path(data = df_ell, aes(x = x, y = y, colour = Genotype), size = 1, linetype = 2) +
    scale_colour_manual(values = pal) +
    labs(title = paste0(rsid, ', All individuals'),
         x = 'Allele A (C) Intensity',
         y = 'Allele B (G) Intensity') +
    theme_bw() +
    theme(legend.key = element_blank())
  
  print(p)
  ggsave(paste0(rsid, '.pdf'), p, width = 10, height = 8)
  
}

# examples
plot.intensity(rsid = rsid[1])
plot.intensity(rsid = rsid[2])
plot.intensity(rsid = rsid[3])
plot.intensity(rsid = rsid[4])
