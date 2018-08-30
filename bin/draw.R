###############plot stripe before #############
rm(list = ls())
setwd("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/plot_data")
stripeData <- read.table("./aggregate_data/aggregate_data.txt", header = T, sep = "\t")
stripeData$C <- as.character(stripeData$C)
colnames(stripeData) <- c("X", "Y", "Stripe")

library(latex2exp)
library(ggplot2)

g <- ggplot(data = stripeData)
g <- g + geom_point(aes( x = X, y = Y, color=Stripe), show.legend = T, alpha=0.8, size=1)

g <- g + theme_linedraw()
g <- g + xlab("GC content") + ylab(TeX('$\\log (D^T/D^N)$')) 
#theme(strip.text = element_text(colour = 'red', size=10.5), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white") )

g <- g + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key = element_rect(fill = "white"),
  legend.title=element_text(size=10.5),
  #legend.position = c(0.5, 0.6),
  #legend.position = "bottom",
  legend.background = element_blank(),
  legend.box.background = element_rect(colour = "black"),
  #legend.direction="horizontal",
  axis.text=element_text(size=10.5),
  axis.title=element_text(size=10.5),
  axis.text.x = element_text(size=10.5),
  axis.text.y = element_text(size=10.5),
  plot.title = element_text(size=10.5))

png(paste("/home/dustin/github/mythesis/figures/stripeAggregate.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()

png(paste("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/figures/stripeAggregate.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()

g
###############Draw BAF#################


rm(list = ls())
setwd("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/plot_data")

BAFClusterCenter <- read.table("./BAF_plot_data/BAFPlot_centers_data.txt", header = T, sep = "\t")
colnames(BAFClusterCenter) <- c("BAF", "Stripe")
BAFClusterCenter$Stripe <- as.character(BAFClusterCenter$Stripe)

BAFData <- read.table("./BAF_plot_data/BAFPlot_data.txt", header = T, sep = "\t")
colnames(BAFData) <- c("BAF", "Class", "Stripe")

BAFData$Stripe <- as.character(BAFData$Stripe)
BAFData$Class <- as.character(BAFData$Class)

library(latex2exp)
library(ggplot2)

g <- ggplot(data = BAFData, aes( x = Stripe, y = BAF));
g <- g + geom_violin()
g <- g + geom_point(data=BAFClusterCenter, aes(x=Stripe, y=BAF), size=1)
g <- g + theme_linedraw()
g <- g + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text=element_text(size=10.5),
  axis.title=element_text(size=10.5),
  axis.text.x = element_text(size=10.5),
  axis.text.y = element_text(size=10.5),
  plot.title = element_text(size=10.5))

png(paste("/home/dustin/github/mythesis/figures/stripeBAFClustering.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()
png(paste("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/figures/stripeBAFClustering.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()
g

################draw after#################
rm(list = ls())
setwd("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/plot_data")
stripeData <- read.table("./decompose_plot_data/decompose_data.txt", header = T, sep = "\t")
stripeData$C <- as.character(stripeData$C)
colnames(stripeData) <- c("X", "Y", "Stripe")

library(latex2exp)
library(ggplot2)

g <- ggplot(data = stripeData)
g <- g + geom_point(aes( x = X, y = Y, color=Stripe), show.legend = F, alpha=0.8, size=1)

g <- g + theme_linedraw()
g <- g + xlab("GC content") + ylab(TeX('$\\log (D^T/D^N)$')) 
#theme(strip.text = element_text(colour = 'red', size=10.5), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white") )

g <- g + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key = element_rect(fill = "white"),
  legend.title=element_text(size=10.5),
  #legend.position = c(0.5, 0.6),
  #legend.position = "bottom",
  legend.background = element_blank(),
  legend.box.background = element_rect(colour = "black"),
  #legend.direction="horizontal",
  axis.text=element_text(size=10.5),
  axis.title=element_text(size=10.5),
  axis.text.x = element_text(size=10.5),
  axis.text.y = element_text(size=10.5),
  plot.title = element_text(size=10.5))

png(paste("/home/dustin/github/mythesis/figures/stripeDecompose.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()

png(paste("/media/d/docEdit/tempArticle/BIBM2017/p-SCNAClonal/figures/stripeDecompose.png", sep = ""), pointsize = 10.5, family = "Times", width = 2.8, height= 2.5, units="in", res=600)
g
dev.off()

g
