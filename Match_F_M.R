## Find the overlap of SNAs between males and females in EMMAX (continuous or binary)
## Date: Feb 22, 2015

snpf = read.csv('Markers LS F.csv', header = TRUE, as.is = TRUE)  # load snp table of female (continuous outcome)
snpm = read.csv('Markers LS M.csv', header = TRUE, as.is = TRUE)  # load snp table of male (continuous outcome)
snpfb = read.csv('Markers LS F_binary.csv', header = TRUE, as.is = TRUE)  # load snp table of female (binarized outcome)
snpmb = read.csv('Markers LS M_binary.csv', header = TRUE, as.is = TRUE)  # load snp table of male (binarized outcome)
###### SNPs selected from lifespan GWAS in EMMAX model
###### Outcome lifespan was coded as a continuous or binary variable

set1 = snpf$Marker  # extract set 1 from table: female snp id (character vector)
set2 = snpm$Marker  # extract set 2 from table: male snp id (character vector)
set3 = snpfb$Marker  # extract set 3 from table: female (binary) snp id (character vector)
set4 = snpmb$Marker  # extract set 4 from table: male (binary) snp id (character vector)

area1 = length(set1)  # get size of set 1
area2 = length(set2)  # get size of set 2
area3 = length(set3)  # get size of set 3
area4 = length(set4)  # get size of set 4

n12 = length(intersect(set1, set2))  # get intersect set size of set 1 and set 2
n13 = length(intersect(set1, set3))  # get intersect set size of set 1 and set 3
n14 = length(intersect(set1, set4))  # get intersect set size of set 1 and set 4
n23 = length(intersect(set2, set3))  # get intersect set size of set 2 and set 3
n24 = length(intersect(set2, set4))  # get intersect set size of set 2 and set 4
n34 = length(intersect(set3, set4))  # get intersect set size of set 3 and set 4
n123 = length(intersect(intersect(set1, set2), set3))  # get intersect set size of set 1, 2 and 3
n124 = length(intersect(intersect(set1, set2), set4))  # get intersect set size of set 1, 2 and 4
n134 = length(intersect(intersect(set1, set3), set4))  # get intersect set size of set 1, 3 and 4
n234 = length(intersect(intersect(set2, set3), set4))  # get intersect set size of set 2, 3 and 4
n1234 = length(intersect(intersect(intersect(set1, set2), set3), set4))  # get intersect set size of set 1, 2, 3 and 4

library(VennDiagram)  # library for drawing venn diagrams in R
# another usable R package is: venneuler
# For set visualization techniques, see the Nat. Met. paper:
# http://www.nature.com/nmeth/journal/v11/n8/full/nmeth.3033.html

# draw venn plot (see the package help for more details)
# we have four sets, so use draw.quad.venn() here
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
  category = c("Female", "Male", "Female Binary", "Male Binary"),  # specify the set names
  fill = c("orange", "red", "green", "blue"),  # specify set colors
  lty = "dashed",  # specify dash pattern of the circles' circumferences
  cex = 2,  # specify the size of the set labels
  cat.cex = 2,  # specify the size of the set names
  cat.col = c("orange", "red", "green", "blue")  # specify the color of the set names
)

png("quad_venn.png", width = 800, height = 600)  # open a png device, output a 800 x 600 png file
grid.draw(venn.plot)  # VennDiagram uses the grid device to plot, so here is a bit different with base graphics. This syntax is directly from the example of the function's help.
dev.off()  # close the device (and the png file will be written to disk)
