## Find the overlap of SNAs between males and females in EMMAX (continuous or binary)
## Date: Feb 22, 2015

snpf = read.csv('Markers LS F.csv', header = TRUE, as.is = TRUE)
snpm = read.csv('Markers LS M.csv', header = TRUE, as.is = TRUE)
snpfb = read.csv('Markers LS F_binary.csv', header = TRUE, as.is = TRUE)
snpmb = read.csv('Markers LS M_binary.csv', header = TRUE, as.is = TRUE)
###### SNPs selected from lifespan GWAS in EMMAX model
###### Outcome lifespan was coded as a continuous or binary variable

# set 1: F
# set 2: M
# set 3: F Binary
# set 4: M Binary

set1 = snpf$Marker
set2 = snpm$Marker
set3 = snpfb$Marker
set4 = snpmb$Marker

area1 = length(set1)
area2 = length(set2)
area3 = length(set3)
area4 = length(set4)

n12 = length(intersect(set1, set2))
n13 = length(intersect(set1, set3))
n14 = length(intersect(set1, set4))
n23 = length(intersect(set2, set3))
n24 = length(intersect(set2, set4))
n34 = length(intersect(set3, set4))
n123 = length(intersect(intersect(set1, set2), set3))
n124 = length(intersect(intersect(set1, set2), set4))
n134 = length(intersect(intersect(set1, set3), set4))
n234 = length(intersect(intersect(set2, set3), set4))
n1234 = length(intersect(intersect(intersect(set1, set2), set3), set4))

library(VennDiagram)

venn.plot <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = c("Female", "Male", "Female Binary", "Male Binary"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)

png("quad_venn.png", width = 800, height = 600)
grid.draw(venn.plot)
dev.off()
