library(ggplot2)
library(gridExtra)
library(sybilSBML)
library(lattice)
library(gridExtra)
library(gtable)
library(grid)
library(gridExtra)


#------------------------------------------------------------------------------#
#      				Ecoli data                          						                       #
#------------------------------------------------------------------------------#


pathEcoli<-"/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/without_constraint/result/Free_met_at_full-knockout/"
#iJO model


print(load(paste0(pathEcoli,"0_at_WT_SUM_Perc_Ecoli_without_const20161101_1750.Rdata")))
plWTE<-EplSum
plWTE[56]<-0# s0001


print(load(paste0(pathEcoli,"0_ATP.Rdata")))
plATPE<-EplSumFree
plATPE[56]<-0

print(load(paste0(pathEcoli,"0_GTP.Rdata")))
plGTPE<-EplSumFree
plGTPE[56]<-0

print(load(paste0(pathEcoli,"0_CTP.Rdata")))
plCTPE<-EplSumFree
plCTPE[56]<-0

print(load(paste0(pathEcoli,"0_ITP.Rdata")))
plITPE<-EplSumFree
plITPE[56]<-0

print(load(paste0(pathEcoli,"0_UTP.Rdata")))
plUTPE<-EplSumFree
plUTPE[56]<-0

print(load(paste0(pathEcoli,"0_NADH.Rdata")))
plNADHE<-EplSumFree
plNADHE[56]<-0

print(load(paste0(pathEcoli,"0_NADPH.Rdata")))
plNADPHE<-EplSumFree
plNADPHE[56]<-0

print(load(paste0(pathEcoli,"0_FADH2.Rdata")))
plFADH2E<-EplSumFree
plFADH2E[56]<-0

print(load(paste0(pathEcoli,"0_Q8H2.Rdata")))
plQ8H2E<-EplSumFree
plQ8H2E[56]<-0

print(load(paste0(pathEcoli,"0_GLU.Rdata")))
plGLUE<-EplSumFree
plGLUE[56]<-0

print(load(paste0(pathEcoli,"0_ACCOA.Rdata")))
plACCOAE<-EplSumFree
plACCOAE[56]<-0

print(load(paste0(pathEcoli,"0_MQL8.Rdata")))
plDMMQL8E<-EplSumFree
plDMMQL8E[56]<-0

print(load(paste0(pathEcoli,"0_MQL8.Rdata")))
plMQL8E<-EplSumFree
plMQL8E[56]<-0



#------------------------------------------------------------------------------#
#      				yeast data                          					   #
#------------------------------------------------------------------------------#


pathyeast<-"/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/without_constraint/results/Free_met_at_full-knockout/"


print(load(paste0(pathyeast,"0_at_WT_SUM.Rdata")))
plWTY<-EplSum

print(load(paste0(pathyeast,"0_ATP.Rdata")))
plATPY<-EplSumFree

print(load(paste0(pathyeast,"0_CTP.Rdata")))
plCTPY<-EplSumFree

print(load(paste0(pathyeast,"0_GTP.Rdata")))
plGTPY<-EplSumFree

print(load(paste0(pathyeast,"0_UTP.Rdata")))
plUTPY<-EplSumFree

print(load(paste0(pathyeast,"0_NADH.Rdata")))
plNADHY<-EplSumFree

print(load(paste0(pathyeast,"0_NADPH.Rdata")))
plNADPHY<-EplSumFree

print(load(paste0(pathyeast,"0_ACCOA.Rdata")))
plACCOAY<-EplSumFree

print(load(paste0(pathyeast,"0_GLU.Rdata")))
plGLUY<-EplSumFree






#------------------------------------------------------------------------------#
#      				FigureS1                           					   #
#------------------------------------------------------------------------------#


Fig1<-function(WT,X,spx){
    
    orrM<-WT!=0 | X!=0
    plWT<-as.numeric(names(table(WT[orrM])))
    countWT<-unname(table(WT[orrM]))
    
    plX<-as.numeric(names(table(X[orrM])))
    countX<-unname(table(X[orrM]))
    
    datsWT <- data.frame(pl =plWT,count= as.integer(countWT),lines="A",stringsAsFactors = FALSE)
    datsX <- data.frame(pl =plX,count= as.integer(countX),lines="B",stringsAsFactors = FALSE)
    
    
    
    
    
    datsWTX<-data.frame(rbind(datsWT,datsX))
    
    maxx<-max(datsWT$pl,datsX$pl)
    maxy<-max(datsWT$count,datsX$count)
    PNE<-ggplot(datsWTX, aes(x = pl, y = count, fill = lines)) +
    geom_bar(stat = "identity", position=position_dodge(),width=0.7)+
    labs(size= "Nitrogen",
    x = " Pleiotropy",
    y = "Count")+ theme_bw()+
    scale_x_continuous(breaks = c(seq(0,maxx,2),maxx))+
    scale_y_continuous(breaks =c(seq(0,maxy,5),maxy))+ theme(legend.position=c(0.95,0.95),legend.justification=c(1,1))+
    ggtitle(spx)+ theme(plot.title = element_text(face="italic",size=20))+
    theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=14),axis.text.y = element_text(face="bold", color="black",
    size=12))+theme(legend.title=element_blank(),axis.title=element_text(size=18),legend.text=element_text(size=18))
    
    PNE<-PNE+ scale_fill_discrete(breaks=c("A", "B"),labels=c("Wildtype", "Free NADPH"))
    
    return(PNE)
}

X<-plNADPHE
WT<-plWTE

PNE<-Fig1(WT,X,"E. coli")

X<-plNADPHY
WT<-plWTY
PNY<-Fig1(WT,X,"S. cerevisiae")



pdf(file="S1.pdf",
width=14,height=10)
grid.arrange(PNE,PNY ,nrow=2)
dev.off()







#------------------------------------------------------------------------------#
#      				FigureS8                           					   #
#------------------------------------------------------------------------------#
#Ecoli
# these results were calculated beasd on
#  a  length(which(WT-FREE_met!=0))

dat10<-data.frame(mets=c("ATP","CTP","GTP","UTP","ITP","NADH","NADPH","FADH2","FMNH2","Q8H2","MQL8","DMMQL8","ACCOA","GLU")
,poc=c(5.14,5.61,5.61,5.61,5.61,6.07,17.29,0,0,1.87,1.87,1.87,6.07,6.07))


PPE<-ggplot(dat10, aes(x=mets,y=poc))+geom_bar(position="dodge", stat="identity",width=.5,fill="#F8766D")+ theme_bw()+
scale_size_area() +
labs(size= "Nitrogen",
x = "Currency metabolites",
y = "Percent of change")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle("E. coli")+ theme(
plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold"))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())+scale_y_continuous(breaks =seq(0,18,2))



#------------------------------------------------------------------------------#
#   #S. cerevisiae   				                            					   #
#------------------------------------------------------------------------------#
# these results were calculated beasd on
#  a  length(which(WT-FREE_met!=0))
dat11<-data.frame(mets=c("ATP","CTP","GTP","UTP","ITP","NADH","NADPH","FADH2","FMNH2","ACCOA","GLU"),
pocy=c(52.40,7.42,15.28,30.57,0,46.29 ,43.67 ,44.99,44.99,27.1,43.67))

PPY<-ggplot(dat11, aes(x=mets,y=pocy))+geom_bar(position="dodge", stat="identity",width=.5,fill="#F8766D")+ theme_bw()+
scale_size_area() +
labs(size= "Nitrogen",
x = "Currency metabolites",
y = "Percent of change")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle("S. cerevisiae")+ theme(
plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold"))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())+scale_y_continuous(breaks =seq(0,52,4))




pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure5/Fig_5_Ecoli_yeast_withoutCoinst.pdf")
grid.arrange(PPE,PPY ,nrow=2)
dev.off()






















##########################################################################

print(load("/Users/deya/Desktop/mount_HPC/HILBERT/iJO/results_without_const/allWT/PLWT_Ecoli_oc.Rdata"))
plWTE<-EplSum
print(load("/Users/deya/Desktop/mount_HPC/HILBERT/iJO/results_without_const/free/PL_NADPH_Ecoli_ow.Rdata"))
plNADPHE<-EplSumFree

for(i in 1:201){
    
    plWTE[[i]][56]<-0
    plNADPHE[[i]][56]<-0
}


print(load("/Users/deya/Desktop/mount_HPC/HILBERT/Yeast_7.6/without_const_results/allWT/PLWT_yeast_oc.Rdata"))
plWTY<-EplSum
print(load("/Users/deya/Desktop/mount_HPC/HILBERT/Yeast_7.6/without_const_results/free/PL_NADPH_yeast_oc.Rdata"))
plNADPHY<-EplSumFree



##########################################################################
###################Figure2 #################################################

WTs<-plWTY
NADPH<-plNADPHY

dd<-cbind(plWT[[1]]-plNADPH[[1]])
GL<-which(dd!=0)
for(i in 1:length(GL)){
    p1<-do.call( rbind, WTs)[,GL[i]]
    p2<-do.call( rbind, NADPH)[,GL[i]]
    
    
    dat8 <- data.frame(dens = c(p1,p2), lines = rep(c("Wildtype", "FREE NADPH"),each=length(p1)),order=rep(c(1,2),each=length(p1)),Value = seq(0,1,0.005))
    
    p<- ggplot(dat8, aes(x=Value, y=dens)) +
    geom_step(size=1.5,aes(group=order, color=lines)) +
    scale_y_continuous(breaks = seq(0,max(dat8[,1]),2))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
    labs(
    x = "Wildtype flux ratio",
    y = "Pleiotropy")+theme_bw()+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=8, angle=45))
    
    fname = paste0("/Users/deya/Desktop/RM_plot/test_yeast_without/",format(Sys.time(), paste0(GL[i],".pdf")))
    
    pdf(file = fname)
    print(p)
    dev.off()
}


WTs<-plWTE
NADPH<-plNADPHE

i=155#iJO
#i=300#iJO
p1<-do.call( rbind, WTs)[,i]
p2<-do.call( rbind, NADPH)[,i]


dat8 <- data.frame(dens = c(p1,p2), lines = rep(c("Wildtype", "FREE NADPH"),each=length(p1)),order=rep(c(1,2),each=length(p1)),Value = seq(0,1,0.005))

dat8$lines<- factor(dat8$lines, levels = rev(levels(dat8$lines)))# reverse order

p<- ggplot(dat8, aes(x=Value, y=dens , order = -as.numeric(lines))) +theme_bw()+
geom_step(size=1.5,aes(group=order, color=lines)) +
scale_y_continuous(breaks = c(seq(0,max(dat8[,1]),4),max(dat8[,1])))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
labs(
x = "Wildtype flux ratio",
y = "Pleiotropy")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=9))+  ggtitle("E. coli")+ theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+theme(legend.title=element_blank())






i=373
#i=188
WTs<-plWTY
NADPH<-plNADPHY


p1<-do.call( rbind, WTs)[,i]
p2<-do.call( rbind, NADPH)[,i]


dat9 <- data.frame(dens = c(p1,p2), lines = rep(c("Wildtype", "FREE NADPH"),each=length(p1)),order=rep(c(1,2),each=length(p1)),Value = seq(0,1,0.005))

dat9$lines<- factor(dat9$lines, levels = rev(levels(dat9$lines)))# reverse order

py<- ggplot(dat9, aes(x=Value, y=dens , order = -as.numeric(lines))) +theme_bw()+
geom_step(size=1.5,aes(group=order, color=lines)) +
scale_y_continuous(breaks = c(seq(0,max(dat9[,1]),2),max(dat9[,1])))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
labs(
x = "Wildtype flux ratio",
y = "Pleiotropy")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=9))+  ggtitle("S. cerevisiae")+ theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+theme(legend.title=element_blank())




pdf(file="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure2/Fig2_without_const_gene_LPD.pdf")
grid.arrange(p,py ,nrow=2)
dev.off()



##########################################################################
###################Figure3 ###############################################

S0S1<-function(plSum){
    
    unq1<-list()
    fitness_Y<-list()
    for(i in 1:length((plSum[[1]]))){
        
        # print(i)
        fitness_Y[[i]]<- do.call( rbind, plSum)[,i]
        
        unq<-unique(fitness_Y[[i]])
        
        unq1[i]<-length(unq)-1
        
        
    }
    
    notM<-unlist(unq1)[which(unlist(unq1)!=0)]
    stepsMM<-unq1
    
    dd<-lapply(fitness_Y[ which(plSum[[1]]!=0)],function(x) which(x[1]==x[200]))
    #length(unlist(dd))
    dd1<-lapply(dd,function(x) is.integer(x) && length(x) == 0L)
    Ste0E<-which(plSum[[1]]!=0)[which(unlist(dd1)==F)]
    #include #Constant pleiotropy
    Stp1C<-which(unlist(unq1)==1)
    Stp1<-length(Stp1C[which(Stp1C %in% Ste0E ==F)])
    Stp0<-length(Stp1C[which(Stp1C %in% Ste0E ==T)])
    
    count<-unname(table(notM))
    steps<-as.numeric(names(table(notM)))
    #browser()
    sq1<-seq(min(steps),max(steps))
    cs1<-which(is.na(match(sq1,steps)))#-1
    if(length(cs1)==0){
        count<-count
        steps<-steps
    }else{
        count<-c(count,rep(0,length(cs1)))
        steps<-c(steps,cs1)
    }
    
    steps<-c(0,steps) #Zero step and pl>0
    count <-c(Stp0,count)
    
    steps[which(steps==1)]<-1
    count[which(steps==1)]<-Stp1
    
    
    
    return(list(count,steps))
}

Fig3<-function(datsMM,datsRM,spx){
    
    
    datz<-rbind(datsRM,datsMM)
    maxx<-max(datsMM$steps,datsRM$steps)
    maxy<-max(datsMM$counts,datsRM$counts)
    PNE<-ggplot(datz, aes(x = steps, y = counts, fill = order)) +
    geom_bar(stat = "identity", position=position_dodge(),width=.6)+
    labs(size= "Nitrogen",
    x = "Number of steps",
    y = "Count")+ theme_bw()+
    scale_x_continuous(breaks = seq(0,maxx))+
    scale_y_continuous(breaks =c(seq(0,maxy,10),maxy))+ theme(legend.position=c(1,1),legend.justification=c(1,1))+  ggtitle(spx)+ theme(plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=10),axis.text.y = element_text(face="bold", color="black",
    size=9))+theme(legend.title=element_blank())
    
    PNE<-PNE+ scale_fill_discrete(breaks=c("A", "B"),labels=c("Wildtype", "Free NADPH"))
    
    return(PNE)
}




TabelE<-S0S1(plWTE)
datsMM <- data.frame(counts =as.integer(TabelE[[1]]),steps=as.integer(TabelE[[2]]),order="A",stringsAsFactors = FALSE)
TabelEN<-S0S1(plNADPHE)
datsRM <- data.frame(counts =as.integer(TabelEN[[1]]),steps=as.integer(TabelEN[[2]]),order="B",stringsAsFactors = FALSE)

PNE<-Fig3(datsMM,datsRM,"E. coli")

TabelY<-S0S1(plWTY)
datsMM <- data.frame(counts =as.integer(TabelY[[1]]),steps=as.integer(TabelY[[2]]),order="A",stringsAsFactors = FALSE)
TabelYN<-S0S1(plNADPHY)
datsRM <- data.frame(counts =as.integer(TabelYN[[1]]),steps=as.integer(TabelYN[[2]]),order="B",stringsAsFactors = FALSE)
PNY<-Fig3(datsMM,datsRM,"S. cerevisiae")
##########Final Fig 3 New##############


#pdf(file="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure3/Fig3_without_const.pdf")
pdf(file="/Users/deya/Dropbox/Paper_Figures_WO_s0001/Fig3_without_const.pdf")
grid.arrange(PNE,PNY ,nrow=2)
dev.off()


##########Final Fig 3 New##############





maxx<-max(datsMM$steps,datsRM$steps)
maxy<-max(datsMM$counts,datsRM$counts)
PNY<-ggplot() +
geom_step(data=datsMM, aes(x=steps, y=counts,colour=order)) +
geom_step(data=datsRM, aes(x=steps, y=counts,colour=order))+
labs(size= "Nitrogen",
x = "Number of steps",
y = "Count")+ theme_bw()+
scale_x_continuous(breaks = seq(0,maxx))+
scale_y_continuous(breaks =seq(0,maxy,10))+ theme(legend.position=c(1,1),legend.justification=c(1,1))+  ggtitle("E. coli")+ theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())

PNE<-PNE+ scale_colour_discrete(breaks=c("A", "B"),labels=c("Wildtype", "Free NADPH"))

PNY<-PNY+ scale_colour_discrete(breaks=c("A", "B"),labels=c("Wildtype", "Free NADPH"))

pdf(file="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure3/Fig3_without_const.pdf")
grid.arrange(PNE,PNY ,nrow=2)
dev.off()


#S. cerevisiae
##########################################################################
###################Figure 4 ###########################################


plsteps<-function(plsums){
    unq1<-list()
    fitness_Y<-list()
    for(i in 1:length((plsums[[1]]))){
        
        fitness_Y[[i]]<- do.call( rbind, plsums)[,i]
        
        unq<-unique(fitness_Y[[i]])
        unq1[[i]]<-length(unq)-1
    }
    stepsMM<-unq1
    
    
    datsMM <- data.frame(pl =plsums[[1]],steps=unlist(stepsMM))
    
    
    dd<-lapply(fitness_Y[ which(plsums[[1]]!=0)],function(x) which(x[1]==x[200]))
    #length(unlist(dd))
    dd1<-lapply(dd,function(x) is.integer(x) && length(x) == 0L)
    Ste0E<-which(plsums[[1]]!=0)[which(unlist(dd1)==F)]
    #include #Constant pleiotropy
    Stp1C<-which(unlist(unq1)==1)
    Stp1<-Stp1C[which(Stp1C %in% Ste0E ==F)]
    Stp0<-Stp1C[which(Stp1C %in% Ste0E ==T)]
    
    
    datsMM[Stp0,2]<-0#Constant pleiotropy
    datsMM[Stp1,2]<-1
    
    return(datsMM)
}

datsMM<-plsteps(plWTE)
datsMM<-plsteps(plWTY)

datsMM<-datsMM[datsMM[,1]!=0,]
maxx<-max(datsMM$pl)
maxy<-max(datsMM$steps)

datsMM1<-count_(datsMM, vars = c("pl","steps"))
freq1<-sort(unlist(unique(datsMM1[,3])))

pMMY <- ggplot(datsMM, aes(pl, steps)) + geom_point()+ theme_bw()+scale_size_area()+
labs(size= "Nitrogen",
x = "pleiotropy at full gene knockout",
y = "Number of steps")+ ggtitle("S. cerevisiae \nWildtype") + theme(plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+
scale_x_continuous(breaks = c(seq(0,maxx,2),maxx))+scale_y_continuous(breaks =c(seq(0,maxy,2),maxy))+
theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())+ stat_sum()+scale_size(breaks=freq1,labels=freq1)

#cor(datsMM[,1],datsMM[,2],method="spearman")


#pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure4/Ecoli_yeast_without_const2.pdf")
pdf(file ="/Users/deya/Dropbox/Paper_Figures_WO_s0001/Ecoli_yeast_without_const2.pdf")
grid.arrange(pMM,pMMY ,nrow=2)
dev.off()





datsRM<-plsteps(plNADPHE)
#datsRM<-plsteps(plNADPHY)

datsRM<-datsRM[datsRM[,1]!=0,]
datsRM<-data.frame(pl=datsRM[,1],steps=datsRM[,2])

datsRM1<-count_(datsRM, vars = c("pl","steps"))

freq1<-sort(unlist(unique(datsRM1[,3])))
maxx<-max(datsRM1$pl)
maxy<-max(datsRM1$steps)


pRMY <- ggplot(datsRM1, aes(pl, steps))+ geom_point(aes(size = n))+ theme_bw()+scale_size_area()+
labs(size= "Nitrogen",
x = "pleiotropy at full gene knockout",
y = "Number of steps")+ ggtitle("S. cerevisiae  \nFree NADPH") + theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+
scale_x_continuous(breaks = c(seq(0,maxx,2),maxx))+scale_y_continuous(breaks =c(seq(0,maxy,2),maxy))+
theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank(),legend.key.size = unit(.4, "cm"))+scale_size(breaks=freq1,labels=freq1)




pdf(file ="/Users/deya/Dropbox/Paper_Figures_WO_s0001/Ecoli_yeast_NADPH_without_const2.pdf")
#pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure4/Ecoli_yeast_NADPH_without_const2.pdf")
grid.arrange(pRM,pRMY ,nrow=2)
dev.off()

#cor(datsRM[,1][datsRM[,1]!=0],datsRM[,2][datsRM[,1]!=0],method="spearman")




#######################################Figure SM ###############################################################



plWT<-plWTE
plWT<-plWTY


plWT<-plNADPHE#Ecoli

plWT<-plNADPHY#Yeast

fitness_Y<-list()
for(i in 1:length(plWT[[1]])){
    
    fitness_Y[[i]]<- do.call( rbind, plWT)[,i]
    
}


unq<-sapply(fitness_Y,function(x) unique(x))
unq1<-sapply(unq,function(x) (length(x)-1))

#dat81 <- data.frame(unlist(fitness_Y), lines = rep(unq1,each=201),Value = seq(0,1,0.005))

DFs<-lapply(names(table(unq1)), function(x) which(x==unq1))
DFs1<-lapply(DFs,function(x) x[1])
dat81 <- data.frame(unlist(fitness_Y[unlist(DFs1)]), lines = rep(unq1[unlist(DFs1)],each=201),Value = seq(0,1,0.005))



colnames(dat81)[1]<-"pl"
colnames(dat81)[2]<-"steps"
stps<-sort(unique(unq1))

p1<- ggplot(dat81, aes(x=Value,y=pl)) +theme_bw()+
geom_step(size=1.2,aes(color=steps,group=steps)) +
scale_y_continuous(breaks = c(seq(0,max(dat81[,1]),4),max(dat81[,1])))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
labs(
x = "Wildtype flux ratio",
y = "Pleiotropy")+ theme(axis.text.x = element_text(face="bold", color="black",
size=8, angle=45,hjust=1))+ ggtitle("E. coli") + theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+ theme(legend.text=element_text(size=10))+ guides(size = guide_legend(order = 3))

pMM<-p1+ scale_colour_gradientn(colours=rainbow(10))

pMY<-p1+ scale_colour_gradientn(colours=rainbow(10))


pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/S2_NADPH_allGenes_Groupsteps_Ecoli_yeast_With_cons.pdf")

grid.arrange(pMM,pMY ,nrow=2)
dev.off()

#S. cerevisiae
#+ theme(legend.position=c(1,1),legend.justification=c(1,1))
#+ theme(legend.key.size=unit(0.5,"cm"))











