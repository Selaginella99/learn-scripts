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
