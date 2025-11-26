#do a subtree from plant megatree with selected species
#he megatree is ‘GBOTB.extended.WP.tre’ reported by Jin and Qian (2022).
devtools::install_github("jinyizju/U.PhyloMaker")
library("U.PhyloMaker")
setwd("~/Documents/Labo/Projets/2019_EC2Co_Orchamp/Primers_STR/Sym_genes_genoms")
sp.list <- read.csv("2_taxa_tree.csv")
megatree <- read.tree('https://raw.githubusercontent.com/megatrees/plant_20221117/main/plant_megatree.tre')
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/plant_20221117/main/plant_genus_list.csv', sep=",")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result
#tree as newick to open in itol
write.tree(result$phylo, "output_tree.tre")
#write table with plant families
write.table(result$sp.list, file="sp_list_fam.txt", quote=FALSE, sep="\t")


