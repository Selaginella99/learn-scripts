## Find the overlap of SNAs of lifespan and changing of the weight between males and females in EMMAX (quantitative outcome)
## Date: March 2, 2015

snpf = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/Markers LS F.csv', header = TRUE, as.is = TRUE)  # load snp of lifespan table of female (continuous outcome)
snpm = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/Markers LS M.csv', header = TRUE, as.is = TRUE)  # load snp of lifespan table of male (continuous outcome)
snpf_w = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/Markers RS F.csv', header = TRUE, as.is = TRUE)  # load snp table of female (continuous outcome)
snpm_w = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/Markers RS M.csv', header = TRUE, as.is = TRUE)  # load snp table of male (continuous outcome)
###### SNPs selected from lifespan GWAS in EMMAX model

set1 = snpf$Marker  # extract set 1 from table: female snp id (character vector)
set2 = snpm$Marker  # extract set 2 from table: male snp id (character vector)
set3 = snpf_w$Marker  # extract set 3 from table: female snp id (character vector)
set4 = snpm_w$Marker  # extract set 4 from table: male snp id (character vector)

area1 = length(set1)  # get size of set 1
area2 = length(set2)  # get size of set 2
area3 = length(set3)  # get size of set 3
area4 = length(set4)  # get size of set 4

n12_lr = length(intersect(set1, set2))  # get intersect set size of set 1 and set 2
n13_lr = length(intersect(set1, set3))  # get intersect set size of set 1 and set 3
n14_lr = length(intersect(set1, set4))  # get intersect set size of set 1 and set 4
n23_lr = length(intersect(set2, set3))  # get intersect set size of set 2 and set 3
n24_lr = length(intersect(set2, set4))  # get intersect set size of set 2 and set 4
n34_lr = length(intersect(set3, set4))  # get intersect set size of set 3 and set 4
n123_lr = length(intersect(intersect(set1, set2), set3))  # get intersect set size of set 1, 2 and 3
n124_lr = length(intersect(intersect(set1, set2), set4))  # get intersect set size of set 1, 2 and 4
n134_lr = length(intersect(intersect(set1, set3), set4))  # get intersect set size of set 1, 3 and 4
n234_lr = length(intersect(intersect(set2, set3), set4))  # get intersect set size of set 2, 3 and 4
n1234_lr = length(intersect(intersect(intersect(set1, set2), set3), set4))  # get intersect set size of set 1, 2, 3 and 4

###
nn12_lr = intersect(set1, set2)  # get intersect set of set 1 and set 2
nn13_lr = intersect(set1, set3)  # get intersect set of set 1 and set 3
nn14_lr = intersect(set1, set4)  # get intersect set of set 1 and set 4
nn23_lr = intersect(set2, set3)  # get intersect set of set 2 and set 3
nn24_lr = intersect(set2, set4)  # get intersect set of set 2 and set 4
nn34_lr = intersect(set3, set4)  # get intersect set of set 3 and set 4
nn123_lr = intersect(intersect(set1, set2), set3)  # get intersect set of set 1, 2 and 3
nn124_lr = intersect(intersect(set1, set2), set4)  # get intersect set of set 1, 2 and 4
nn134_lr = intersect(intersect(set1, set3), set4)  # get intersect set of set 1, 3 and 4
nn234_lr = intersect(intersect(set2, set3), set4)  # get intersect set of set 2, 3 and 4
nn1234_lr = intersect(intersect(intersect(set1, set2), set3), set4)  # get intersect set of set 1, 2, 3 and 4
###

library(VennDiagram)  # library for drawing venn diagrams in R
# another usable R package is: venneuler
# For set visualization techniques, see the Nat. Met. paper:
# http://www.nature.com/nmeth/journal/v11/n8/full/nmeth.3033.html

# draw venn plot (see the package help for more details)
# we have four sets, so use draw.quad.venn() here
venn.plot_lr <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12_lr,
  n13 = n13_lr,
  n14 = n14_lr,
  n23 = n23_lr,
  n24 = n24_lr,
  n34 = n34_lr,
  n123 = n123_lr,
  n124 = n124_lr,
  n134 = n134_lr,
  n234 = n234_lr,
  n1234 = n1234_lr,
  category = c("female lifespan", "male lifespan", "female weight", "male weight"),  # specify the set names
  fill = c("orange", "red", "green", "blue"),  # specify set colors
  lty = "dashed",  # specify dash pattern of the circles' circumferences
  cex = 1.78,  # specify the size of the set labels
  cat.cex = 1.78,  # specify the size of the set names
  cat.col = c("orange", "red", "green", "blue")  # specify the color of the set names
)

png("quad_venn_lr_emmax.png", width = 800, height = 600)  # open a png device, output a 800 x 600 png file
grid.draw(venn.plot_lr)  # VennDiagram uses the grid device to plot, so here is a bit different with base graphics. This syntax is directly from the example of the function's help.
dev.off()  # close the device (and the png file will be written to disk)

write.csv(nn12_lr,file="Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/LS_emmax.csv")
LS.emmax = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/LS_emmax.csv', header = TRUE, as.is = TRUE)
write.csv(nn34_lr,file="Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/RS_binary.csv")
RS.emmax = read.csv('Z:/Project/Tools/Golden Helix/Training/Data_test_FRAM/IBC GRU + NPU/IBC49K/RS_emmax.csv', header = TRUE, as.is = TRUE)

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

n12_ls = intersect(set1, set2)  # get intersect set of set 1 and set 2
n13_ls = intersect(set1, set3)  # get intersect set of set 1 and set 3
n14_ls = intersect(set1, set4)  # get intersect set of set 1 and set 4
n23_ls = intersect(set2, set3)  # get intersect set of set 2 and set 3
n24_ls = intersect(set2, set4)  # get intersect set of set 2 and set 4
n34_ls = intersect(set3, set4)  # get intersect set of set 3 and set 4
n123_ls = intersect(intersect(set1, set2), set3)  # get intersect set of set 1, 2 and 3
n124_ls = intersect(intersect(set1, set2), set4)  # get intersect set of set 1, 2 and 4
n134_ls = intersect(intersect(set1, set3), set4)  # get intersect set of set 1, 3 and 4
n234_ls = intersect(intersect(set2, set3), set4)  # get intersect set of set 2, 3 and 4
n1234_ls = intersect(intersect(intersect(set1, set2), set3), set4)  # get intersect set of set 1, 2, 3 and 4


## Plot the intersect sets of set1, set2, set3 and set 4 using Venn plot
library(VennDiagram)  # library for drawing venn diagrams in R
# another usable R package is: venneuler
# For set visualization techniques, see the Nat. Met. paper:
# http://www.nature.com/nmeth/journal/v11/n8/full/nmeth.3033.html

# draw venn plot (see the package help for more details)
# we have four sets, so use draw.quad.venn() here
venn.plot_lf <- draw.quad.venn(
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
  category = c("Female LS Quant", "Male LS Quant", "Female LS Binary", "Male LS Binary"),  # specify the set names
  fill = c("orange", "red", "green", "blue"),  # specify set colors
  lty = "dashed",  # specify dash pattern of the circles' circumferences
  cex = 1.7,  # specify the size of the set labels
  cat.cex = 1.7,  # specify the size of the set names
  cat.col = c("orange", "red", "green", "blue")  # specify the color of the set names
)

## draw.venn.plot: http://127.0.0.1:13913/library/VennDiagram/html/draw.quad.venn.html

png("quad_venn_ls.png", width = 800, height = 600)  # open a png device, output a 800 x 600 png file
grid.draw(venn.plot_lf)  # VennDiagram uses the grid device to plot, so here is a bit different with base graphics. This syntax is directly from the example of the function's help.
dev.off()  # close the device (and the png file will be written to disk)
