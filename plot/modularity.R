library(ggplot2)
library(gridExtra)
library(sybilSBML)
library(lattice)
library(gridExtra)
library(gtable)
library(grid)
library(gridExtra)
library(reshape2)
require(gplots)
library("RColorBrewer")
library(gridGraphics)
library(grid)
library(cobs)#concave



#iJO model
print(load("data_plot/modularity/plWTE1.Rdata"))

mycols = c("#DEEBF7","#4291C2")



pdf(file="Ecoli.pdf",height=10)
heatmap.2(plWTE1, scale = "none", col = mycols, main = "E. coli",cexRow =0.22,cexCol = 0.7, key=FALSE, trace="none")
dev.off()

#levelplot(plWTE1,scales=list(x=list(at=seq(1,10),labels=labXY),y=list(at=seq(1,11),
#labels=labXY)),border="black",main="Fitenss vs. pl in Min", xlab="Fintenss", ylab="pl",panel = myPanel,col.regions=c("white",rainbow(10)))



#------------------------------------------------------------------------------#
#      				histogram                           					   #
#------------------------------------------------------------------------------#


EQbs = read.table("data_plot/modularity//Qbs.Ecoli.10000.txt",header = T)
#YGqb<-0.1971

Edata <- data.frame(Value=EQbs)


pE<- ggplot(Edata, aes(EQbs) ) + theme_bw()+ geom_histogram(fill="#F8766D",bins=20,binwidth=0.002)+
scale_x_continuous(limits = c(0.10,0.25),breaks = seq(0.10,0.25,0.01))+
scale_y_continuous(breaks = seq(0,1500,100))+
scale_size_area() +labs(size= "Nitrogen",
x = "Modularity",
y = "Count")+ ggtitle("E. coli")+
theme(plot.title = element_text(size=20, face="italic",hjust = 0))+
theme(axis.text.x = element_text(face="bold", color="black",
size=14),axis.text.y = element_text(face="bold", color="black",
size=14),axis.title=element_text(size=16))


pdf(file="Yeast_hist.pdf",width=12,height=12)
grid.arrange(pE)
dev.off()

#0.2348 observed value





#------------------------------------------------------------------------------#
#      				                            					   #
#------------------------------------------------------------------------------#




#yeast model

print(load("plWTY1.RData"))

mycols = c("#DEEBF7","#4291C2")

print(load("plWTE1.Rdata"))
pdf(file="east.pdf",height=10)
heatmap.2(plWTY1, scale = "none", col = mycols, main = "S. cerevisiae",cexRow =0.3,cexCol = 0.8, key=FALSE, trace="none")
dev.off()



#------------------------------------------------------------------------------#
#      				histogram                           					   #
#------------------------------------------------------------------------------#




YQbs = read.table("data_plot/modularity/Qbs.yeast.10000.txt",header = T)
#YGqb<-0.1971

data <- data.frame(Value=YQbs)


pY<- ggplot(data, aes(YQbs) ) + theme_bw()+ geom_histogram(fill="#F8766D",bins=20,binwidth=0.002)+
scale_x_continuous(limits = c(0.06,0.20),breaks = seq(0.06,0.20,0.01))+
scale_y_continuous(breaks = seq(0,1000,100))+
scale_size_area() +labs(size= "Nitrogen",
x = "Modularity",
y = "Count")+ ggtitle("S. cerevisiae")+
theme(plot.title = element_text(size=22, face="italic",hjust = 0))+
theme(axis.text.x = element_text(face="bold", color="black",
size=14),axis.text.y = element_text(face="bold", color="black",
size=14),axis.title=element_text(size=16))


pdf(file="Yeast_hist.pdf",width=12,height=12)
grid.arrange(pE)
dev.off()



#------------------------------------------------------------------------------#
#      				                           					   #
#------------------------------------------------------------------------------#


pdf(file="Ecoli_Yeast_Modularity1_row.pdf",width=16,height=8)
#grid.arrange(pE,pY,nrow=2)
grid.arrange(pE,pY,ncol=2)
dev.off()


#------------------------------------------------------------------------------#
#      				                           					   #
#------------------------------------------------------------------------------#



pdf(file="Ecoli_yeast.pdf",height=10)
for(i in 1:2){
    heatmap.2(plWTE1, scale = "none", col = mycols, main = "E. coli",cexRow =0.22,cexCol = 0.7, key=FALSE, trace="none")
    ,heatmap.2(plWTY1, scale = "none", col = mycols, main = "S. cerevisiae",cexRow =0.33,cexCol = 0.8, key=FALSE, trace="none")
    , ncol=2, nrow =1)
    
    dev.off()
    
    
    
    par(mfrow = c(2, 2))
    heatmap.2(plWTE1, scale = "none", col = mycols, main = "E. coli",cexRow =0.22,cexCol = 0.7, key=FALSE, trace="none")
    heatmap.2(plWTY1, scale = "none", col = mycols, main = "S. cerevisiae",cexRow =0.33,cexCol = 0.8, key=FALSE, trace="none")
    par(mfrow = c(1, 1))
    
    
    
    
    
    
    #------------------------------------------------------------------------------#
    #      				Figure                           					   #
    #------------------------------------------------------------------------------#
    
    library(gridGraphics)
    library(grid)
    
    grab_grob <- function(){
        grid.echo()
        grid.grab()
    }
    
    arr <- replicate(2, matrix(sample(1:100),nrow=10,ncol=10), simplify = FALSE)
    arr[[1]]<-plWTE1
    arr[[1]]<-plWTE1
    spec<-c("E. coli","S. cerevisiae")
    library(gplots)
    gl <- lapply(1:4, function(i){
        heatmap.2(arr[[i]] ,dendrogram ='row',
        Colv=FALSE, col=greenred(800),
        key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
        trace='none', colsep=1:10,
        sepcolor='white', sepwidth=0.05,
        scale="none",cexRow=0.2,cexCol=2,
        labCol = colnames(arr[[i]]),
        hclustfun=function(c){hclust(c, method='mcquitty')},
        lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),
        )
        grab_grob()
    })
    
    
    grid.newpage()
    library(gridExtra)
    grid.arrange(grobs=gl, ncol=2, clip=TRUE)
    
    
    source("heatmap.21.R")
    library(gtools)
    
    
    heatmap.21(plWTE1)
    heatmap.21(plWTY1)
    
    
    
    library(gridGraphics)
    library(grid)
    
    grab_grob <- function(){
        grid.echo()
        grid.grab()
    }
    
    matx<-list(plWTE1,plWTY1)
    
    library(gplots)
    gl <- lapply(1:2, function(i){
        heatmap.2(matx[[i]])
        grab_grob()
    })
    
    grid.newpage()
    library(gridExtra)
    grid.arrange(grobs=gl, ncol=2, clip=TRUE)
    
    
    
    
    library(gridGraphics)
    library(grid)
    heatmap(as.matrix(mtcars))
    
    library(gridGraphics)
    grab_grob <- function(){
        grid.echo()
        grid.grab()
    }
    
    g <- grab_grob()
    grid.newpage()
    
    
    matx<-list(plWTE1,plWTY1)
    
    spec<-c("E. coli","S. cerevisiae")
    library(gplots)
    gl <- lapply(1:2, function(i){
        
        
        heatmap.2(matx[[i]])#, scale = "none", col = mycols, main =spec[i] ,cexRow =0.22,cexCol = 0.7, key=FALSE, trace="none") 
        
        grab_grob()
    })
    
    grid.newpage()
    library(gridExtra)
    grid.arrange(grobs=gl, ncol=2, clip=TRUE)
    
    
    
    
    
    
