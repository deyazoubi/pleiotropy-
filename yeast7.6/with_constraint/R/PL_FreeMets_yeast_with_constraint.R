#free_model<-Emodel

#------------------------------------------------------------------------------#
#           allow free currancy mets. in the model                             #
#------------------------------------------------------------------------------#
cm<-c("ATP","CTP","GTP","UTP","ITP","NADH","NADPH","FADH2","FMNH2","GLU","ACCOA")

currency_met <- function (mod3,indx) {
    
    #cm<- readline("Please enter your Free currency metabolite?")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print(cm[indx])
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    print("#------------------------------------------------------------------------------#")
    
    
    if(cm[indx]=="ATP"){
        
        free_model<-addReact(mod3,id="ATPM",met=met_id(mod3)[c(324,614,288,605,1045)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="CTP"){
        
        free_model<-addReact(mod3,id="CTPM",met=met_id(mod3)[c(410,614,350,605,1045)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="GTP"){
        
        free_model<-addReact(mod3,id="GTPM",met=met_id(mod3)[c(596,614,560,605,1045)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="UTP"){
        
        free_model<-addReact(mod3,id="UTPM",met=met_id(mod3)[c(1214,614,1193,605,1045)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="ITP"){
        
        free_model<-addReact(mod3,id="ITPM",met=met_id(mod3)[c(751,614,649,605,1045)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="NADH"){
        
        free_model<-addReact(mod3,id="NADHM",met=met_id(mod3)[c(957,605,952)],Scoef=c(-1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="NADPH"){
        
        free_model<-addReact(mod3,id="NADPHM",met=met_id(mod3)[c(965,605,961)],Scoef=c(-1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="FADH2"){
        
        free_model<-addReact(mod3,id="FADH2M",met=met_id(mod3)[c(532,605,530)],Scoef=c(-1,2,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="FMNH2"){
        
        free_model<-addReact(mod3,id="FMNH2M",met=met_id(mod3)[c(547,605,544)],Scoef=c(-1,2,1),reversible =T,lb=-1000,ub=1000)
        
    }else if(cm[indx]=="GLU"){
        
        free_model<-addReact(mod3,id="GLUM",met=met_id(mod3)[c(792,614,124,310,605)],Scoef=c(-1,-1,1,1,2),reversible =T,lb=-1000,ub=1000)#L-Glutamate
        
    }else if(cm[indx]=="ACCOA"){
        
        free_model<-addReact(mod3,id="ACCOAM",met=met_id(mod3)[c(614,274,605,266,400)],Scoef=c(-1,-1,1,1,1),reversible =T,lb=-1000,ub=1000)
        
    }else{
        
        stop("missing currency metabolite")
        
    }
    
    EbioFree<-biomass_reactsADDALL1(free_model,2984)
    
    return(EbioFree)
    
}




#------------------------------------------------------------------------------#
#wt production of the orignal model to compare them with components productions(free case)  #
#------------------------------------------------------------------------------#


wtProd <- function (comp,model) {
    
    
    mod<-changeObjFunc(model, react=comp,obj_coef=1)
    wt1<-optimizeProb(mod)
    wt<-lp_obj(invisable_values(wt1 ))
    
    return(wt)
    
    
}

wts <- mapply ( wtProd,RNX_id, MoreArgs =list(Ebio[[1]]), SIMPLIFY = TRUE )




#------------------------------------------------------------------------------#
#                              calculating pleiotropy                          #
#------------------------------------------------------------------------------#

# The following function simulates the effect of all possible gene deletions
# in the production of particular biomass components


affectedBiomassComponent <- function (wt,comp,model,j,imp) {
    
    
    
    mod<-changeObjFunc(model, react=comp,obj_coef=1)
    sgds <- onegenedel(mod,alg="fba1",geneList =allGenes(mod)[imp],flux_Percent=j)
    ogd<-lp_obj(invisable_values(sgds ))
    
    red <- (abs(wt)-abs(ogd))/abs(wt) >= 1E-4
    
    return ( red )
}


run_pl<-function(wt,comp,model,fluxpercent,imp){
    
    # pl <- mapply ( affectedBiomassComponent,wt,comp, MoreArgs =list(model,fluxpercent,imp ), SIMPLIFY = TRUE )
    pl <- mcmapply ( affectedBiomassComponent,wt,comp, MoreArgs =list(model,fluxpercent,imp ), SIMPLIFY = TRUE,mc.cores=2)
    return(pl)
}


# currency metabolites loop
for(k in 1: 11) {
    
    FreeMod<-currency_met(Emodel,k)
    
    ReacExF<-FreeMod[[2]]
    
    RNX_idF<-match(ReacExF,react_id(FreeMod[[1]]))
    
    
    
    
    #------------------------------------------------------------------------------#
    #simulated mutations of different severity by restricting the flux through
    #all reactions catalyzed by the gene to a fixed fraction of the wildtype flux,
    #which I reduced in 0.5%-steps
    #------------------------------------------------------------------------------#
    
    for(f in 1:length(fluxPercent)){
        
        EplSumFree<-rowSums(run_pl(unname(wts),RNX_idF,FreeMod[[1]],fluxPercent[f],seq(1,length(allGenes(Emodel)))))
        
        fname = format(Sys.time(), paste0(paste0(f,"_"),paste0(cm[k],".Rdata")))
        save(EplSumFree,file=paste0("results/",fname))
        save.image()
        
        rm(EplSumFree)
        unlink(".RData")
        
    }
    
}




