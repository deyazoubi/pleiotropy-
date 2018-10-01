library(ggplot2)
library(gridExtra)
library(sybilSBML)
library(lattice)
library(gridExtra)
library(gtable)
library(grid)
library(gridExtra)
library(cobs)#concave


#------------------------------------------------------------------------------#
#      				Ecoli data                          						                       #
#------------------------------------------------------------------------------#


pathEcoli<-"/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/with_constraint/result/Free_met_at_full-knockout/"
#iJO model
print(load(paste0(pathEcoli,"0_at_WT_SUM_Perc_Ecoli20161101_1512.Rdata")))
plWTE<-EplSum

print(load(paste0(pathEcoli,"0_ATP.Rdata")))
plATPE<-EplSumFree

print(load(paste0(pathEcoli,"0_CTP.Rdata")))
plCTPE<-EplSumFree

print(load(paste0(pathEcoli,"0_GTP.Rdata")))
#new
plGTPE<-EplSumFree

print(load(paste0(pathEcoli,"0_UTP.Rdata")))
#new
plUTPE<-EplSumFree

print(load(paste0(pathEcoli,"0_NADH.Rdata")))
plNADHE<-EplSumFree

print(load(paste0(pathEcoli,"0_NADPH.Rdata")))
plNADPHE<-EplSumFree

print(load(paste0(pathEcoli,"0_FADH2.Rdata")))
plFADH2E<-EplSumFree

print(load(paste0(pathEcoli,"0_Q8H2.Rdata")))
plQ8H2E<-EplSumFree

print(load(paste0(pathEcoli,"0_ACCOA.Rdata")))
#new
plACCOAE<-EplSumFree

print(load(paste0(pathEcoli,"0_GLU.Rdata")))
plGLUE<-EplSumFree





#------------------------------------------------------------------------------#
#      				yeast data                          					   #
#------------------------------------------------------------------------------#


pathyeast<-"/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/with_constraint/result/Free_met_at_full-knockout/"


#yeast model
print(load(paste0(pathyeast,"0_at_WT_yeast20161104_1518.Rdata")))

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
#      				Figure1                           					   #
#------------------------------------------------------------------------------#

#match("s0001",allGenes(Ecoli))
#56

Fig1<-function(datsWT,datsX,spx){
    
    orrM<-WT!=0 | X!=0
    #exclude "s0001" gene
    if(length(WT)>1300){
        
        orrM[56]<-F
        
    }
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
    scale_x_continuous(breaks = seq(0,maxx,1))+
    scale_y_continuous(breaks =seq(0,maxy,4))+ theme(legend.position=c(0.95,0.95),legend.justification=c(1,1))+  ggtitle(spx)+
    theme(plot.title = element_text(face="italic",size=20))+
    theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=14, angle=45),axis.text.y = element_text(face="bold", color="black",
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


pdf(file="Figure1_NADPH_140317.pdf",width=10,hight=5)
grid.arrange(PNE,PNY ,nrow=2)
dev.off()





#########################Figure2 and 3  S6 S7 data#################################################


print(load("/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/with_constraint/result/PLWT_Ecoli.Rdata"))
plWTE<-EplSum

print(load("/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/with_constraint/result/PL_NADPH_Ecoli.Rdata"))
plNADPHE<-EplSumFree



print(load("/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/with_constraint/result/PLWT_yeast.Rdata"))
plWTY<-EplSum

print(load("/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/with_constraint/result/PL_NADPH_yeast.Rdata"))
plNADPHY<-EplSumFree



#------------------------------------------------------------------------------#
#      				Figure2                           					   #
#------------------------------------------------------------------------------#

i=155#iJO

WTs<-plWTE
NADPH<-plNADPHE

p1<-do.call( rbind, WTs)[,i]
p2<-do.call( rbind, NADPH)[,i]


dat8 <- data.frame(dens = c(p1,p2), lines = rep(c("Wildtype", "FREE NADPH"),each=length(p1)),order=rep(c(1,2),each=length(p1)),Value = seq(0,1,0.005))

dat8$lines<- factor(dat8$lines, levels = rev(levels(dat8$lines)))# reverse order


p<- ggplot(dat8, aes(x=Value, y=dens, order = -as.numeric(lines))) +theme_bw()+
geom_step(size=1.5,aes(group=order, color=lines)) +
scale_y_continuous(breaks = c(seq(0,max(dat8[,1]),4),max(dat8[,1])))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
labs(
x = "Wildtype flux ratio",
y = "Pleiotropy")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=9))+  ggtitle("E. coli")+ theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+theme(legend.title=element_blank())


x<-dat8[,4]
y<-dat8[,1]
conreg(x=x,y=y)



i=373#yeast

WTs<-plWTY
NADPH<-plNADPHY


p1<-do.call( rbind, WTs)[,i]
p2<-do.call( rbind, NADPH)[,i]


dat9 <- data.frame(dens = c(p1,p2), lines = rep(c("Wildtype", "FREE NADPH"),each=length(p1)),order=rep(c(1,2),each=length(p1)),Value = seq(0,1,0.005))

dat9$lines<- factor(dat9$lines, levels = rev(levels(dat9$lines)))# reverse order


py<- ggplot(dat9, aes(x=Value, y=dens, order = -as.numeric(lines))) +theme_bw()+
geom_step(size=1.5,aes(group=order, color=lines)) +
scale_y_continuous(breaks = c(seq(0,max(dat9[,1]),4),max(dat9[,1])))+ scale_x_continuous(breaks = seq(0,1,0.05))+scale_size_area()+
labs(
x = "Wildtype flux ratio",
y = "Pleiotropy")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=9))+  ggtitle("S. cerevisiae")+ theme(plot.title = element_text(face="italic"))+
theme(plot.title = element_text(hjust = 0))+theme(legend.title=element_blank())



pdf(file="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure2/Fig2_WC_gene_ATP1.pdf")
grid.arrange(p,py ,nrow=2)
dev.off()






#------------------------------------------------------------------------------#
#      				Figure3                           					   #
#------------------------------------------------------------------------------#

# histogram of the # of steps include zero steps MM####

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
    #browser()
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
    scale_y_continuous(breaks =seq(0,80,4))+ theme(legend.position=c(0.95,0.95),legend.justification=c(1,1))+  ggtitle(spx)+ theme(plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=10),axis.text.y = element_text(face="bold", color="black",
    size=10))+theme(legend.title=element_blank())
    
    PNE<-PNE+ scale_fill_discrete(breaks=c("A", "B"),labels=c("Wildtype", "Free NADPH"))
    
    return(PNE)
}


for(i in 1:201){
    
    plWTE[[i]][56]<-0
    plNADPHE[[i]][56]<-0
}

#do.call( rbind, plWTE)[,56]
#do.call( rbind, plNADPHE)[,56]


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






pdf(file="/Users/deya/Dropbox/Paper_Figures_WO_s0001/Fig3_WC.pdf")
grid.arrange(PNE,PNY ,nrow=2)
dev.off()






#------------------------------------------------------------------------------#
#      				Figure4                           					   #
#------------------------------------------------------------------------------#
#Ecoli
dat10<-data.frame(mets=c("ATP","CTP","GTP","UTP","ITP","NADH","NADPH","FADH2","FMNH2","Q8H2","MQL8","DMMQL8","ACCOA","GLU"),
# these results were calculated beasd on
#  a  length(which(WT-FREE_met!=0))
poc=c(8.7,4.7,5.7, 4.7, 0 ,41.3,10.7,31.7,0,39.0,0,0,4.7,30.7))

PPE<-ggplot(dat10, aes(x=mets,y=poc))+geom_bar(position="dodge", stat="identity",width=.5,fill="#F8766D")+ theme_bw()+
scale_size_area() +
labs(size= "Nitrogen",
x = "Currency metabolites",
y = "Percent of change")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle("E. coli")+ theme(
plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold"))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())+scale_y_continuous(breaks =c(seq(0,max(dat10$poc),5),floor(max(dat10$poc))))
########################################################################################


#S. cerevisiae
dat11<-data.frame(mets=c("ATP","CTP","GTP","UTP","ITP","NADH","NADPH","FADH2","FMNH2","ACCOA","GLU"),
# these results were calculated beasd on
#  a  length(which(WT-FREE_met!=0))
pocy=c(53.60,35.14,13.5,21.6,0.00 ,12.6 ,31.5 ,0,0,1,31.5))

PPY<-ggplot(dat11, aes(x=mets,y=pocy))+geom_bar(position="dodge", stat="identity",width=.5,fill="#F8766D")+ theme_bw()+
scale_size_area() +
labs(size= "Nitrogen",
x = "Currency metabolites",
y = "Percent of change")+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle("S. cerevisiae")+ theme(
plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold"))+ theme(axis.text.x = element_text(face="bold", color="black",
size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
size=10))+theme(legend.title=element_blank())+scale_y_continuous(breaks =c(seq(0,max(dat11$pocy),10),floor(max(dat11$pocy))))


pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure5/Fig_5_Ecoli_yeast_02112016.pdf")
grid.arrange(PPE,PPY ,nrow=2)
dev.off()






#------------------------------------------------------------------------------#
#      				Figure  S2    Ecoli                       					   #
#------------------------------------------------------------------------------#


CMs<-list(WT= plWTE,"Free ATP"=plATPE  ,"Free CTP"=plCTPE,
"Free GTP"=plGTPE, "Free UTP"=plUTPE ,"Free NADH"=plNADHE,
"Free NADPH"= plNADPHE,"Free FADH2"=plFADH2E,"Free Q8H2"=plQ8H2E ,
"Free ACCOA"=plACCOAE,"Free GLU"=plGLUE)


PP<-list()
for(i in 1:(length(CMs)-1)){
    
    
    orrM<-CMs[[1]]!=0 | CMs[[i+1]]!=0
    
    if(length(CMs[[1]])>1300){
        
        orrM[56]<-F
        
    }
    
    plWT<-as.numeric(names(table(CMs[[1]][orrM])))
    countWT<-unname(table(CMs[[1]][orrM]))
    
    plX<-as.numeric(names(table(CMs[[i+1]][orrM])))
    countX<-unname(table(CMs[[i+1]][orrM]))
    
    datsWT <- data.frame(pl =plWT,count= as.integer(countWT),lines="A",stringsAsFactors = FALSE)
    datsX <- data.frame(pl =plX,count= as.integer(countX),lines="B",stringsAsFactors = FALSE)
    datsWTX<-data.frame(rbind(datsWT,datsX))
    
    My<-max(countWT,countX)
    Mx<-max(plWT,plX)
    
    if(i==1){
        ggt<-"E. coli"
    }else{
        ggt<-""
    }
    
    
    if(i==(length(CMs)-1)){
        Xa<-"Pleiotropy"
    }else{
        Xa<-""
    }
    
    
    if(i==4){
        Ya<-"Count"
    }else{
        Ya<-""
    }
    
    PP[[i]]<-ggplot(datsWTX, aes(x = pl, y = count, fill = lines)) +
    geom_bar(stat = "identity", position=position_dodge(),width=0.7)+theme_bw()+
    scale_x_continuous(breaks = c(seq(0,Mx,2),Mx))+
    scale_y_continuous(breaks = seq(0,My,10))+
    scale_size_area() +
    labs(size= "Nitrogen",
    x = Xa,
    y = Ya)+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle(ggt)+ theme(
    plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold",size=10))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=10),axis.text.y = element_text(face="bold", color="black",
    size=10))+labs(fill="")+scale_fill_discrete(breaks=c("A", "B"),labels=c("Wildtype", names(CMs[i+1])))
}


pdf(file="Desktop/S2_Ecoli.pdf",width=7.5,height=30)

grid.arrange(arrangeGrob(PP[[1]],PP[[2]],PP[[3]],PP[[4]],PP[[5]],PP[[7]],PP[[8]],PP[[9]],PP[[10]], ncol=1),heights=unit(1, "npc"))

dev.off()


#------------------------------------------------------------------------------#
#      				Figure  S2    yeast                       					   #
#------------------------------------------------------------------------------#
CMs<-list(WT= plWTY,"Free ATP"=plATPY  , "Free GTP"=plGTPY , "Free CTP"=plCTPY, "Free UTP"=plUTPY,"Free NADH"=plNADHY,"Free NADPH"= plNADPHY, "Free ACCOA"=plACCOAY ,"Free GLU"=plGLUY)

PPY<-list()
for(i in 1:(length(CMs)-1)){
    
    
    orrM<-CMs[[1]]!=0 | CMs[[i+1]]!=0
    
    
    plWT<-as.numeric(names(table(CMs[[1]][orrM])))
    countWT<-unname(table(CMs[[1]][orrM]))
    
    plX<-as.numeric(names(table(CMs[[i+1]][orrM])))
    countX<-unname(table(CMs[[i+1]][orrM]))
    
    datsWT <- data.frame(pl =plWT,count= as.integer(countWT),lines="A",stringsAsFactors = FALSE)
    datsX <- data.frame(pl =plX,count= as.integer(countX),lines="B",stringsAsFactors = FALSE)
    datsWTX<-data.frame(rbind(datsWT,datsX))
    
    My<-max(countWT,countX)
    Mx<-max(plWT,plX)
    
    
    
    
    if(i==1){
        ggt<-"S. cerevisiae"
    }else{
        ggt<-""
    }
    
    
    if(i==8){
        Xa<-"Pleiotropy"
        print(i)
    }else{
        Xa<-""
    }
    
    
    
    if(i==(3)){
        Ya<-"Count"
        print(i)
    }else{
        Ya<-""
    }
    
    PPY[[i]]<-ggplot(datsWTX, aes(x = pl, y = count, fill = lines)) +
    geom_bar(stat = "identity", position=position_dodge(),width=0.7)+theme_bw()+
    scale_x_continuous(breaks = c(seq(1,Mx,2),Mx))+
    scale_y_continuous(breaks = seq(0,My,8))+
    scale_size_area() +
    labs(size= "Nitrogen",
    x = Xa,
    y = Ya)+ theme(legend.position=c(1,1),legend.justification=c(1,1))+ ggtitle(ggt)+ theme(
    plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+theme(legend.position=c(0.90,1),legend.justification=c(1,1),legend.key.size = unit(.4, "cm"),legend.text = element_text(face="bold",size=10))+ theme(axis.text.x = element_text(face="bold", color="black",
    size=10),axis.text.y = element_text(face="bold", color="black",
    size=10))+labs(fill="")+scale_fill_discrete(breaks=c("A", "B"),labels=c("Wildtype", names(CMs[i+1])))
    
}


pdf(file="Desktop/S2_yeast.pdf",width=7.5,height=20)

grid.arrange(arrangeGrob(PPY[[1]],PPY[[2]],PPY[[3]],PPY[[4]],PPY[[5]],PPY[[7]],PPY[[8]], ncol=1),heights=unit(1, "npc"))

dev.off()





#------------------------------------------------------------------------------#
#      				Figure  S6
#please load data line 164			                       					   #
#------------------------------------------------------------------------------#




plotS6<-function(ORG){
    
    fitness_Y<-list()
    for(i in 1:length(plWT[[1]])){
        
        fitness_Y[[i]]<- do.call( rbind, plWT)[,i]
        
    }
    
    
    unq<-sapply(fitness_Y,function(x) unique(x))
    unq1<-sapply(unq,function(x) (length(x)-1))
    
    
    DFs<-lapply(names(table(unq1)), function(x) which(x==unq1))
    DFs1<-lapply(DFs,function(x) x[1])#random
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
    size=8, angle=45,hjust=1))+ ggtitle(ORG) + theme(plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+ theme(legend.text=element_text(size=10))+ guides(size = guide_legend(order = 3))
    
    pMM<-p1+ scale_colour_gradientn(colours=rainbow(10))
    
    return(pMM)
    #pMY<-p1+ scale_colour_gradientn(colours=rainbow(10))
    
    
}


# reload files in lines 164
plWT<-plWTE#Ecoli
plWT<-plNADPHE#Ecoli
pMM<-plotS6("E. coli")


plWT<-plWTY#Yeast
plWT<-plNADPHY#Yeast

pMY<-plotS6("S. cerevisiae")


pdf(file ="Groupsteps_Ecoli_yeast_WC.pdf")
grid.arrange(pMM,pMY ,nrow=2)
dev.off()

#S. cerevisiae
#+ theme(legend.position=c(1,1),legend.justification=c(1,1))
#+ theme(legend.key.size=unit(0.5,"cm"))





library(dplyr)

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



#------------------------------------------------------------------------------#
#      				Figure  S7    
#please load data line 164			                       					   #
#------------------------------------------------------------------------------#


plotS7<-function(datsMM,ORG){
    
    datsMM<-datsMM[datsMM[,1]!=0,]
    datsMM<-data.frame(pl=datsMM[,1],steps=datsMM[,2])
    
    datsMM1<-count_(datsMM, vars = c("pl","steps"))
    
    freq1<-sort(unlist(unique(datsMM1[,3])))
    maxx<-max(datsMM1$pl)
    maxy<-max(datsMM1$steps)
    
    pMM <- ggplot(datsMM1, aes(pl, steps)) + geom_point(aes(size = n))+ theme_bw()+scale_size_area()+
    labs(size= "Nitrogen",
    x = "pleiotropy at full gene knockout",
    y = "Number of steps")+ ggtitle(paste0(ORG,"\nWildtype")) + theme(plot.title = element_text(face="italic"))+theme(plot.title = element_text(hjust = 0))+
    scale_x_continuous(breaks = c(seq(0,maxx,2),maxx))+scale_y_continuous(breaks =c(seq(0,maxy,2),maxy))+
    theme(axis.text.x = element_text(face="bold", color="black",
    size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
    size=10))+theme(legend.title=element_blank(),legend.key.size = unit(.4, "cm"))+scale_size(breaks=freq1,labels=freq1)
    
    
    stat_sum(aes(group = legs))
    
    return(pMM)
    #cor(datsMM[,1][datsMM[,1]!=0],datsMM[,2][datsMM[,1]!=0],method="spearman")
}

datsMM<-plsteps(plWTE)# Ecoli
pMM<-plotS7(datsMM,"E. coli")
datsMMY<-plsteps(plWTY)# yeast
pMMY<-plotS7(datsMMY,"S. cerevisiae")

#pdf(file ="/Users/deya/Dropbox/Pleiotropy&Epistasis/paper/Paper Figures/Figure4/Ecoli_yeast_WT_WC2.pdf")

pdf(file ="/Users/deya/Dropbox/Paper_Figures_WO_s0001/Ecoli_yeast_WT_WC2.pdf")

grid.arrange(pMM,pMMY ,nrow=2)
dev.off()




#------------------------------------------------------------------------------#
#      				Figure  S7   (NADPH)
#please load data line 164			                       					   #
#------------------------------------------------------------------------------#



plotS7_NADPH<-function(datsRM,ORG){
    
    
    datsRM<-datsRM[datsRM[,1]!=0,]
    datsRM<-data.frame(pl=datsRM[,1],steps=datsRM[,2])
    
    datsRM1<-count_(datsRM, vars = c("pl","steps"))
    
    freq1<-sort(unlist(unique(datsRM1[,3])))
    maxx<-max(datsRM1$pl)
    maxy<-max(datsRM1$steps)
    
    
    
    pRM <- ggplot(datsRM1, aes(pl, steps))+ geom_point(aes(size = n))+ theme_bw()+scale_size_area()+
    labs(size= "Nitrogen",
    x = "pleiotropy at full gene knockout",
    y = "Number of steps")+ ggtitle(paste0(ORG, "\nFree NADPH")) + theme(plot.title = element_text(face="italic"))+
    theme(plot.title = element_text(hjust = 0))+
    scale_x_continuous(breaks = c(seq(0,maxx,2),maxx))+scale_y_continuous(breaks =c(seq(0,maxy,2),maxy))+
    theme(axis.text.x = element_text(face="bold", color="black",
    size=10, angle=45, hjust = 1),axis.text.y = element_text(face="bold", color="black",
    size=10))+theme(legend.title=element_blank(),legend.key.size = unit(.4, "cm"))+scale_size(breaks=freq1,labels=freq1)
    
    return(pRM)
    
}


datsRM<-plsteps(plNADPHE)

pRM<-plotS7_NADPH(datsRM,"E. coli")


datsRMY<-plsteps(plNADPHY)
pRMY<-plotS7_NADPH(datsRMY,"S. cerevisiae")



pdf(file ="Ecoli_yeast_NADPH.pdf")
grid.arrange(pRM,pRMY ,nrow=2)
dev.off()

#cor(datsRM[,1][datsRM[,1]!=0],datsRM[,2][datsRM[,1]!=0],method="spearman")



















