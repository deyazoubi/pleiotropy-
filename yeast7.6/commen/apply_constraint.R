

# replace invisable Fitness values with Zero
#=================infeasible Sol. data set <- 0  ====================================
invisable_values<-function(opt){
    
    df <- data.frame(lp_obj(opt), sapply(lp_stat(opt), getMeanStatus))
    infes<-df[df[,2]=="infeasible",]
    
    if(nrow(infes)!=0){
        
        inf_pos<-as.integer(rownames(infes))
        lp_obj(opt)[inf_pos]<-0
    }
    
    
    
    return (opt)
}


#============== apply the constraint on the minmial media  ===========================
constraint_minmial<-function(model){
    
    
    
    fba<-optimizeProb(model,solverParm=list(CPX_PARAM_EPRHS=1e-07))# tolerance to bounds
    mtf <- optimizeProb(model, algorithm = "mtf",wtobj = mod_obj(fba),solverParm=list(CPX_PARAM_EPRHS=1e-07))# tolerance to bounds)))
    wt_flux<-getFluxDist(mtf)
    
    
    model<-constraint_model(model,wt_flux)
    #model<-constraint_model1(model,wt_flux)
    return(model)
}


#============== constraint the model  uppbnd positive and lowbnd negitive wt_flux ===
constraint_model<-function(model,wt_flux){
    #exclude the EX_reaction form the constarint
    ex<-findExchReact(model)
    
    lb<-lowbnd(model)[react_pos(ex)]
    ub<-uppbnd(model)[react_pos(ex)]
    
    ####
    #0<=vi<=viwt   vi>=0
    lowbnd(model)[which( wt_flux >= 0)]<-0
    uppbnd(model)[which( wt_flux >= 0)]<-wt_flux[which( wt_flux >= 0)]
    
     #0>=vj>=vjwt   vi<0
     lowbnd(model)[which( wt_flux < 0)]<-wt_flux[which( wt_flux < 0)]
      uppbnd(model)[which( wt_flux < 0)]<-0
      
    #exclude the EX_reaction form the constarints
    lowbnd(model)[react_pos(ex)]<-lb
    uppbnd(model)[react_pos(ex)]<-ub
    ##
    return(model)
}
