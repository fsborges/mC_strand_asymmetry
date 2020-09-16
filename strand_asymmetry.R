
## Required packages

library(MASS)
library(ggplot2)
library(grid)

sample="wt sc CHG chr1"

## Program

data <- read.table(file="../data/wt_sc_chr1_chg.txt", sep="\t",col.names=c("watson","crick"))

filter <- subset(data,watson>0.25 & crick>0.25)

x <- filter$watson
y <- filter$crick

DF <- data.frame(x,y)

dens <- kde2d(x,y)
gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
names(gr) <- c("xgr", "ygr", "zgr")
mod <- loess(zgr~xgr*ygr, data=gr)
DF$density <- predict(mod, newdata=data.frame(xgr=x, ygr=y))

# Number of filtered cytosines
nrow <- grobTree(textGrob(paste("n =",nrow(filter)),x=0.45,y=0.1,hjust=0,gp=gpar(col="white",fontsize=7,fontface="bold")))
# Average mC with highest density
DF <- DF[order(-DF$density),]
hd <- grobTree(textGrob(head(round(((DF$x+DF$y)/2)*100,1)),x=0.05,y=0.90,hjust=0,gp=gpar(col="white",fontsize=7,fontface="bold")))

# Plot
ggplot(DF, aes(x=x,y=y, color=density)) + 
  geom_point(size=1,shape=16,alpha=1/6) + 
  scale_colour_gradientn(colours=c("darkblue","white","red"),guide = "colourbar") + 
  theme_bw() + 
  labs(x="mC % Watson strand",y="mC % Crick strand",title=sample) +
  theme(panel.background = element_rect(fill = "darkblue"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=8),title =element_text(size=6)) + 
  scale_x_continuous(limits=c(0.25,1.0),breaks = c(0.5,1.0),labels = c("50","100")) + 
  scale_y_continuous(limits=c(0.25,1.0),breaks = c(0.5,1.0),labels = c("50","100")) +
  annotation_custom(nrow) +
  annotation_custom(hd)

# Save
ggsave("../output/wt_sc_CHG_chr1.png",width=2.3,height=1.6,device=function(...)png(...,units='in',res=600))
  