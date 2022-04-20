#Plot trees as distribution

library(phytools)

new_trees <- readRDS("./data/tidy/new_trees_matched_spectra.rds")

#write.tree(new_trees, "./data/tidy/new_trees_matched_spectra.txt")

trees<-readRDS("./data/for_analysis/final_trees_matched_spectra.rds")

class(trees) <- "multiPhylo"

#first attempts
densityTree(trees[1], alpha = 0.25, colors = "blue", nodes = "intermediate", ylim=c(0, length(labels)+1)) #type = "cladogram", nodes = "intermediate",, change font size #, nodes="intermediate",type="chronogram"

plotTree(trees[[1]], colour = col, ylim=c(0, length(labels)+1), type="phylogram", fsize = 0.6)#, tips=tips[trees[[i]]$tip.label], add=FALSE,  )



#try manual way

# tips <- lapply(labels, function(x,y) sapply(y, 
#                                             function(y,x) if(x%in%y$tip.label) which(y$tip.label==x) else NA,
#                                             x=x), y=trees)
# tips <- setNames(sapply(tips, mean, na.rm=TRUE), labels)
# tips <- setNames(1:length(labels), labels[order(tips)])


#make plot 
labels <- trees[[1]]$tip.label
col<-make.transparent("blue",0.01)
jpeg("./output/tree_distribution_final.jpg", height = 10, width = 8, res = 600, units = "in")
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