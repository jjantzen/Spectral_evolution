#Plot trees as distribution

library(phytools)
library(treeio)
library(ggtree)
library(deeptime)

consensus_tree <- readRDS("./data/for_analysis/myc_consensus_tree.rds")

consensus_beast_tree <- read.beast("./data/trees/30burnin_mcc_consensus.tre")

consensus_beast_tree_nex <- read.newick("./data/trees/rooted_mcc_beast_supportvalues.nwk")

try_again <- ape::read.tree("./data/trees/30burnin_mcc_consensus.tre")

trees <- readRDS("./data/for_analysis/myc_tree_92sp_for_analysis.rds")

#class(trees) <- "multiPhylo"

class(beast_trees)

#first attempts
densityTree(trees[1], alpha = 0.25, colors = "blue", nodes = "intermediate", ylim=c(0, length(labels)+1)) #type = "cladogram", nodes = "intermediate",, change font size #, nodes="intermediate",type="chronogram"

plotTree(consensus_tree, colour = col, type="phylogram", fsize = 0.6)#, tips=tips[trees[[i]]$tip.label], add=FALSE,  )

consensus_beast_tree_nex



#try manual way

# tips <- lapply(labels, function(x,y) sapply(y, 
#                                             function(y,x) if(x%in%y$tip.label) which(y$tip.label==x) else NA,
#                                             x=x), y=trees)
# tips <- setNames(sapply(tips, mean, na.rm=TRUE), labels)
# tips <- setNames(1:length(labels), labels[order(tips)])


#make plot 
labels <- trees[[1]]$tip.label
col<-make.transparent("blue",0.01)
jpeg("./output/tree_distribution_myc_trees.jpg", height = 10, width = 8, res = 600, units = "in")
plot.new()
#plot trees - can't change the appearance - how to do that?
for(i in 1:length(trees)){
  if(i==1){
    plotTree(trees[[i]], color = col, ylim=c(0, length(labels)+1), type="phylogram", fsize = 0.6) #tips=tips[trees[[i]]$tip.label],
    par(fg="transparent")
  } else {
    plotTree(trees[[i]], color = col, ylim=c(0, length(labels)+1), add=TRUE, type="phylogram", fsize = 0.6) #tips=tips[trees[[i]]$tip.label], 
  }
}
dev.off()

#set background??
par(fg="black")


max(nodeHeights(trees$STATE_241090000))
#comparison of specific trees (with arrows)
compare.chronograms(new_trees[[1]], new_trees[[2]])
#densityTree(trees,use.edge.length=FALSE,type="phylogram",nodes="inner")


#trim beast consensus and plot

consensus_beast_tree@phylo$tip.label <- gsub("_", " ", consensus_beast_tree@phylo$tip.label)
consensus_beast_tree@phylo$tip.label[which(consensus_beast_tree@phylo$tip.label == "Acer saccharum subsp. nigrum")] <- "Acer nigrum"
consensus_beast_tree@phylo$tip.label[which(consensus_beast_tree@phylo$tip.label == "Alnus incana subsp. rugosa")] <- "Alnus incana"

exclude <- consensus_beast_tree@phylo$tip.label[-which(consensus_beast_tree@phylo$tip.label %in% trees[[1]]$tip.label)]

consensus_beast_tree_pruned <- drop.tip(consensus_beast_tree, tip = exclude)

consensus_beast_tree_pruned@phylo <- ladderize(consensus_beast_tree_pruned@phylo, right = TRUE)

plot(consensus_beast_tree_pruned@phylo, type="phylogram", cex = 0.6, show.tip.label = TRUE, show.node.label = TRUE) #colour = col, 

labels <- consensus_beast_tree_pruned@phylo$tip.label

#plot


p <- ggtree(consensus_beast_tree_pruned)+
  #geom_tiplab(size=3, color="black", aes(label=paste0('italic(', labels, ")")), parse=TRUE)
  geom_tiplab(size=3, color="black", fontface = 3)+
  geom_range(range='height_0.95_HPD', color='red', size=2, alpha=.3)+
  geom_nodelab(aes(x=branch, label=round(posterior, 2), subset=(posterior < 0.996)), vjust=-0.3, size=3, hjust = 0.9)+#
  #xlim(NA, 30)+
  scale_x_continuous(labels = abs(seq(-400,30,100)), limits = c(-400, 30))+#, limits = c(-500, 100))+#, limits = c(-500,100))+
  theme_tree2()#plot.margin = unit(c(14,8,14,8), "mm")
  

revts(p)

jpeg("./output/consensus_tree_ages_support_values.jpg", width = 15, height = 11, units = "in", res = 600)
#pdf("./output/consensus_tree_ages_support_values.pdf", width = 15.5, height = 11) 
revts(p)
dev.off()

#just get rid of negative sign in front of numbers on x axis 



#plotTree(consensus_beast_tree_nex, colour = col, type="phylogram", fsize = 0.6)

str(consensus_beast_tree_pruned)

names(consensus_beast_tree_pruned@data)
