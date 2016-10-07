#17AD 13control 17 PSP

#setwd('../')
library("igraph")
library("rgl")
library("bigmemory")
library("tidyr")
library("corrplot")
library("psych")
nodelocations <- read.csv("all_xys.csv", header=FALSE)
ADminusPSP <- read.csv("ADminusPSP.csv", header=FALSE)
controlminusPSP <- read.csv("controlminusPSP.csv", header=FALSE)
controlminusAD <- read.csv("controlminusAD.csv", header=FALSE)
ConnADPvals <- read.table("ConnADPvals.txt", header=FALSE)
ConnPSPPvals <- read.table("ConnPSPPvals.txt", header=FALSE)
ADPSPPvals <- read.table("ADPSPPvals.txt", header=FALSE)
AvConnMatrix <- read.table("AvConnMatrix.txt", header=FALSE)
AvADMatrix <- read.table("AvADMatrix.txt", header=FALSE)
AvPSPMatrix <- read.table("AvPSPMatrix.txt", header=FALSE)

nodelocations = as.matrix(nodelocations)
ADminusPSP = as.matrix(ADminusPSP)
controlminusPSP = as.matrix(controlminusPSP)
controlminusAD = as.matrix(controlminusAD)
ConnADPvals = as.matrix(ConnADPvals)
ConnPSPPvals = as.matrix(ConnPSPPvals)
ADPSPPvals = as.matrix(ADPSPPvals)
AvConnMatrix = as.matrix(AvConnMatrix)
AvADMatrix = as.matrix(AvADMatrix)
AvPSPMatrix = as.matrix(AvPSPMatrix)

ConnADPvals[is.na(ConnADPvals)] = NA
ConnPSPPvals[is.na(ConnPSPPvals)] = NA
ADPSPPvals[is.na(ADPSPPvals)] = NA
AvConnMatrix[is.na(AvConnMatrix)] = 1
AvADMatrix[is.na(AvADMatrix)] = 1
AvPSPMatrix[is.na(AvPSPMatrix)] = 1

#ConmatrixdifferentfromAD <- cortest.jennrich(AvConnMatrix,AvADMatrix,n1=13, n2=17)  #the Jennrich test for Matrix equality
#ConmatrixdifferentfromPSP <- cortest.jennrich(AvConnMatrix,AvPSPMatrix,n1=13, n2=17)  #the Jennrich test for Matrix equality
#PSPdifferentfromAD <- cortest.jennrich(AvPSPMatrix,AvADMatrix,n1=17, n2=17)  #the Jennrich test for Matrix equality
# 
# ConnADPadjust <- as.matrix(p.adjust(ConnADPvals, method = "fdr"))
# ConnPSPPadjust <- as.matrix(p.adjust(ConnPSPPvals, method = "fdr"))
# ADPSPPadjust <- as.matrix(p.adjust(ADPSPPvals, method = "fdr"))

# corrmatrixtoplot = controlminusPSP Re-order by brain region
regionlist <- read.csv("Atlas_ROIs.csv", header=FALSE)
corrmatrixtoplot <- controlminusPSP[order((regionlist$V2)),order((regionlist$V2))]
occurences<-table(unlist(regionlist$V2))
lobenames <- names(occurences)
clustsizes = NULL

for(i in (1:length(occurences))){
clustsizes[i] = occurences[[lobenames[i]]][1]
}
diag(corrmatrixtoplot) <- 0
#x11()
#corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "Controls minus PSP", tl.pos = "n")
#corrRect(clustsizes,lwd=2)
#dev.copy2pdf(file = "conminusPSPcorrplot.pdf")
#dev.off()   

ConnPSPPvals <- ConnPSPPvals[order((regionlist$V2)),order((regionlist$V2))]
corrmatrixtoplot[ConnPSPPvals >= 0.001] <- 0

plotthesecorrmatrices <- function(x, together=FALSE, group1name = 'Group 1', group2name = 'Group2') {
  # Plots a correlation matrix x then draws two 3d weighted graphs for positive and negative weights if together == FALSE, or combined if together == TRUE
  
  #Optional visualise the connectivity matrix you are plotting (next two lines)
  #x11()
  #corrplot(x, method = 'color', is.corr=FALSE)
  
  # graphtoplot = graph_from_adjacency_matrix(controlminusPSP, weighted = TRUE, diag = FALSE)
  
  graphtoplot <- graph.adjacency(x,weighted=TRUE,mode="lower")
  #E(graphtoplot)[ abs(weight) < max(abs(quantile(controlminusPSP, 0.005, na.rm = TRUE)),abs(quantile(controlminusPSP, 0.995, na.rm = TRUE))) ]$weight <- 0
  graphtoplot <- graphtoplot - E(graphtoplot)[E(graphtoplot)$weight==0]
  
  if(together == TRUE) {
    E(graphtoplot)[ weight > 0]$color <- "blue"
    E(graphtoplot)[ weight < 0]$color <- "red"
    
    rgl.open()
    rglplot(graphtoplot, layout = nodelocations, vertex.size=3, vertex.label=NA, edge.width=3*E(graphtoplot)$weight, edge.color= E(graphtoplot)$color)
    par3d(windowRect = c(20, 30, 800, 800))
    rgl.viewpoint(180,60,90)
    par3d(zoom=0.6)
    legend3d("topright", legend = paste('Both ways for', c(group1name), 'and', c(group2name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
    rgl.snapshot(filename=paste('Both ways for', c(group1name), 'and', c(group2name), '.png'),fmt="png")

  } else {
    tempgraphtoplot <- graphtoplot
    tempgraphtoplot <- delete.edges(tempgraphtoplot, E(tempgraphtoplot)[weight < 0])
    E(tempgraphtoplot)[ weight > 0]$color <- "blue"

    rgl.open()
    rglplot(tempgraphtoplot, layout = nodelocations, vertex.size=3, vertex.label=NA, edge.width=3*E(tempgraphtoplot)$weight, edge.color= E(tempgraphtoplot)$color)
    par3d(windowRect = c(20, 30, 800, 800))
    rgl.viewpoint(180,60,90)
    par3d(zoom=0.6)
    legend3d("topright", legend = paste('Stronger in', c(group1name), 'than', c(group2name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
    rgl.snapshot(filename=paste('Stronger in', c(group1name), 'than', c(group2name), '.png'),fmt="png")
    
    tempgraphtoplot <- graphtoplot
    tempgraphtoplot <- delete.edges(tempgraphtoplot, E(tempgraphtoplot)[weight > 0])
    E(tempgraphtoplot)[ weight < 0]$color <- "red"

    rgl.open()
    rglplot(tempgraphtoplot, layout = nodelocations, vertex.size=3, vertex.label=NA, edge.width=3*E(tempgraphtoplot)$weight, edge.color= E(tempgraphtoplot)$color)
    par3d(windowRect = c(20, 30, 800, 800))
    rgl.viewpoint(180,60,90)
    par3d(zoom=0.6)
    legend3d("topright", legend = paste('Stronger in', c(group2name), 'than', c(group1name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
    rgl.snapshot(filename=paste('Stronger in', c(group2name), 'than', c(group1name), '.png'),fmt="png")
  }
}

plotthesecorrmatrices(corrmatrixtoplot, group1name = 'Controls', group2name = 'PSP')

# Now repeat for control minus AD
corrmatrixtoplot <- controlminusAD[order((regionlist$V2)),order((regionlist$V2))]
occurences<-table(unlist(regionlist$V2))
lobenames <- names(occurences)
clustsizes = NULL
for(i in (1:length(occurences))){
  clustsizes[i] = occurences[[lobenames[i]]][1]
}
diag(corrmatrixtoplot) <- 0
#x11()
#corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "Controls minus AD", tl.pos = "n")
#corrRect(clustsizes,lwd=2)
#dev.copy2pdf(file = "conminusADcorrplot.pdf")
#dev.off()   

ConnADPvals <- ConnADPvals[order((regionlist$V2)),order((regionlist$V2))]
corrmatrixtoplot[ConnADPvals >= 0.001] <- 0
plotthesecorrmatrices(corrmatrixtoplot, group1name = 'Controls', group2name = 'AD')

#Now repeat for AD minus PSP
corrmatrixtoplot <- ADminusPSP[order((regionlist$V2)),order((regionlist$V2))]
occurences<-table(unlist(regionlist$V2))
lobenames <- names(occurences)
clustsizes = NULL
for(i in (1:length(occurences))){
  clustsizes[i] = occurences[[lobenames[i]]][1]
}
diag(corrmatrixtoplot) <- 0
#x11()
#corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "AD minus PSP", tl.pos = "n")
#corrRect(clustsizes,lwd=2)
#dev.copy2pdf(file = "ADminusPSPcorrplot.pdf")
#dev.off()   

ADPSPPvals <- ADPSPPvals[order((regionlist$V2)),order((regionlist$V2))]
corrmatrixtoplot[ADPSPPvals >= 0.001] <- 0
plotthesecorrmatrices(corrmatrixtoplot, group1name = 'AD', group2name = 'PSP')

# Now make simpler graphs based on unified brain regions

plotthesecontractedcorrmatrices <- function(corrmatrixtoplot, group1name = 'Group 1', group2name = 'Group2') {
  diag(corrmatrixtoplot) <- 0
  graphtoplot <- graph.adjacency(corrmatrixtoplot,weighted=TRUE,mode="lower")
  regionlist$V3 <- NA
  
  num = 0
  for(region in unique(unlist(regionlist$V1))){
    num = num+1 
    regionlist$V3[regionlist$V1 == region] <- num
  }
  graphtoplot <- set.vertex.attribute(graphtoplot, 'zlocation', value=nodelocations[,3])
  graphtoplot <- set.vertex.attribute(graphtoplot, 'ylocation', value=nodelocations[,2])
  graphtoplot <- set.vertex.attribute(graphtoplot, 'xlocation', value=nodelocations[,1])
  
  contractedgraphtoplot <- contract.vertices(graphtoplot, regionlist$V3, vertex.attr.comb=list("weight"="mean", "xlocation"="mean", "ylocation"="mean", "zlocation"="mean", "ignore"))
  newlocations <- matrix(c(get.vertex.attribute(contractedgraphtoplot, "xlocation"), get.vertex.attribute(contractedgraphtoplot, "ylocation"), get.vertex.attribute(contractedgraphtoplot, "zlocation")), nrow = 120, ncol = 3)
  E(contractedgraphtoplot)[ abs(weight) < max(abs(quantile(get.edge.attribute(contractedgraphtoplot, "weight"), 0.005, na.rm = TRUE)),abs(quantile(get.edge.attribute(contractedgraphtoplot, "weight"), 0.995, na.rm = TRUE))) ]$weight <- 0
  contractedgraphtoplot <- contractedgraphtoplot - E(contractedgraphtoplot)[E(contractedgraphtoplot)$weight==0]
  
  E(contractedgraphtoplot)[ weight > 0]$color <- "blue"
  E(contractedgraphtoplot)[ weight < 0]$color <- "red"
  
  rgl.open()
  rglplot(contractedgraphtoplot, layout = newlocations, vertex.size=3, vertex.label=NA, edge.width=E(contractedgraphtoplot)$weight, edge.color= E(contractedgraphtoplot)$color)
  par3d(windowRect = c(20, 30, 800, 800))
  rgl.viewpoint(180,60,90)
  par3d(zoom=0.6)
  legend3d("topright", legend = paste('Both ways for', c(group1name), 'and', c(group2name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
  
  tempgraphtoplot <- contractedgraphtoplot
  tempgraphtoplot <- delete.edges(tempgraphtoplot, E(tempgraphtoplot)[weight < 0])
  E(tempgraphtoplot)[ weight > 0]$color <- "blue"
  
  rgl.open()
  rglplot(tempgraphtoplot, layout = newlocations, vertex.size=3, vertex.label=NA, edge.width=3*E(tempgraphtoplot)$weight, edge.color= E(tempgraphtoplot)$color)
  par3d(windowRect = c(20, 30, 800, 800))
  rgl.viewpoint(180,60,90)
  par3d(zoom=0.6)
  legend3d("topright", legend = paste('Stronger in', c(group1name), 'than', c(group2name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
  
  tempgraphtoplot <- contractedgraphtoplot
  tempgraphtoplot <- delete.edges(tempgraphtoplot, E(tempgraphtoplot)[weight > 0])
  E(tempgraphtoplot)[ weight < 0]$color <- "red"
  
  rgl.open()
  rglplot(tempgraphtoplot, layout = newlocations, vertex.size=3, vertex.label=NA, edge.width=3*E(tempgraphtoplot)$weight, edge.color= E(tempgraphtoplot)$color)
  par3d(windowRect = c(20, 30, 800, 800))
  rgl.viewpoint(180,60,90)
  par3d(zoom=0.6)
  legend3d("topright", legend = paste('Stronger in', c(group2name), 'than', c(group1name)), pch = 16, col = c("red"), cex=4, inset=c(0.02))
  
}

corrmatrixtoplot <- controlminusPSP[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'Control', group2name = 'PSP')

corrmatrixtoplot <- controlminusAD[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'Control', group2name = 'ADMCI')

corrmatrixtoplot <- ADminusPSP[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'ADMCI', group2name = 'PSP')

corrmatrixtoplot <- AvConnMatrix[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'Control Pos', group2name = 'Control Neg')
x11()
corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "Controls", tl.pos = "n")
corrRect(clustsizes,lwd=2)
dev.copy2pdf(file = "RawAvConPlot.pdf")
dev.off()  

corrmatrixtoplot <- AvADMatrix[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'AD Pos', group2name = 'AD Neg')
x11()
corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "AD", tl.pos = "n")
corrRect(clustsizes,lwd=2)
dev.copy2pdf(file = "RawAvADPlot.pdf")
dev.off()  

corrmatrixtoplot <- AvPSPMatrix[order((regionlist$V2)),order((regionlist$V2))]
plotthesecontractedcorrmatrices(corrmatrixtoplot, group1name = 'PSP Pos', group2name = 'PSP Neg')
x11()
corrplot(corrmatrixtoplot, method = 'color', is.corr=FALSE, title = "PSP", tl.pos = "n")
corrRect(clustsizes,lwd=2)
dev.copy2pdf(file = "RawAvPSPPlot.pdf")
dev.off()  


# 
# 
#
# Old code that I really should delete now:
#
# diag(controlminusAD) <- 0
# x11()
# corrplot(controlminusAD, method = 'color', is.corr=FALSE)
# graphtoplot <- graph.adjacency(controlminusAD,weighted=TRUE,mode="lower")
# 
# # graphtoplot = graph_from_adjacency_matrix(controlminusAD, weighted = TRUE, diag = FALSE)
# 
# E(graphtoplot)[ abs(weight) < max(abs(quantile(controlminusAD, 0.005, na.rm = TRUE)),abs(quantile(controlminusAD, 0.995, na.rm = TRUE))) ]$weight <- 0
# graphtoplot <- graphtoplot - E(graphtoplot)[E(graphtoplot)$weight==0]
# E(graphtoplot)[ weight < 0]$color <- "blue"
# E(graphtoplot)[ weight > 0]$color <- "red"
# 
# rgl.open()
# rglplot(graphtoplot, layout = nodelocations, vertex.size=2, vertex.label=NA, edge.width=E(graphtoplot)$weight, edge.color= E(graphtoplot)$color)
# legend3d("topright", legend = paste('Stronger in', c('Control', 'AD')), pch = 16, col = c("red", "blue"), cex=1, inset=c(0.02))
# 
# diag(ADminusPSP) <- 0
# x11()
# corrplot(ADminusPSP, method = 'color', is.corr=FALSE)
# graphtoplot <- graph.adjacency(ADminusPSP,weighted=TRUE,mode="lower")
# 
# # graphtoplot = graph_from_adjacency_matrix(ADminusPSP, weighted = TRUE, diag = FALSE)
# 
# E(graphtoplot)[ abs(weight) < max(abs(quantile(ADminusPSP, 0.005, na.rm = TRUE)),abs(quantile(ADminusPSP, 0.995, na.rm = TRUE))) ]$weight <- 0
# graphtoplot <- graphtoplot - E(graphtoplot)[E(graphtoplot)$weight==0]
# E(graphtoplot)[ weight < 0]$color <- "blue"
# E(graphtoplot)[ weight > 0]$color <- "red"
# 
# rgl.open()
# rglplot(graphtoplot, layout = nodelocations, vertex.size=2, vertex.label=NA, edge.width=E(graphtoplot)$weight, edge.color= E(graphtoplot)$color)
# legend3d("topright", legend = paste('Stronger in', c('AD', 'PSP')), pch = 16, col = c("red", "blue"), cex=1, inset=c(0.02))


#### Read in single subject data and do stats on that to work out significantly changed connection strengths
## AD - subject 001-017, control 018-030, PSP 031-047
# 
# filenames = list.files(path="../DatatosendtoTim/data_for_GraphVar/matrices/", pattern="results+.*txt")
# varnames <- strsplit(filenames,'.txt')
# subject_numbers <- extract_numeric(filenames)
# 
# for(i in 1:47){
#   filepath <- file.path("../DatatosendtoTim/data_for_GraphVar/matrices",paste(varnames[i],'.txt',sep=""))
#   #assign(i, read.table(filepath, header=FALSE))
#   if(subject_numbers[i]<18000){ ## UP TO HERE - MAKE AN ARRAY BASED ON GROUP
#     AllADconnmatrices <- lapply(filepath, read.table, header=FALSE)
#     do.call(rbind, AllADconnmatrices)
#   }
#   
#     
#     
# }
# 
# 


