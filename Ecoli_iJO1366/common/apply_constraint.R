constraint_model1<-function(model,wt_flux){
    browser()
    # alternative:
    
    NaffeReact<-which(gpr(model)=="")
    
    lb<-lowbnd(model)[NaffeReact]
    ub<-uppbnd(model)[NaffeReact]
    
    mtf1.current <- abs(wt_flux)
    uppbnd(model) <- mtf1.current
    rev<-react_rev(model)
    idx <- which(rev == TRUE)
    #lowbnd(model)[idx] <- (mtf1.current[idx])*-1
    
    lowbnd(model)[idx] <-ifelse(lowbnd(model)[idx] < ((mtf1.current[idx])*-1),((mtf1.current[idx])*-1),lowbnd(model)[idx])
    
    lowbnd(model)[NaffeReact]<-lb
    uppbnd(model)[NaffeReact]<-ub
    
    return(model)
}


#================ model constraint =================================================
#  use this constraint to use the wt fluxes lb and ub for influx reactions( down regulated flux)
constraint_model_reg<-function(model,wt_flux){
    
    # alternative:
    mtf1.current <- abs(wt_flux)
    uppbnd(model) <- mtf1.current
    rev<-react_rev(model)
    idx <- which(rev == TRUE)
    lowbnd(model)[idx] <- (mtf1.current[idx])*-1
    
    #lowbnd(model)[which(rev == F)]<-0
    return(model)
}

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


#============== constraint the model  uppbnd positive and lowbnd negitive wt_flux ===
constraint_model<-function(model,wt_flux){
    #exclude the EX_reaction form the constarints
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

#============== apply the constraint on the minmial media  ===========================
constraint_minmial<-function(model){
    
    
    
    fba<-optimizeProb(model,solverParm=list(CPX_PARAM_EPRHS=1e-07))# tolerance to bounds
    mtf <- optimizeProb(model, algorithm = "mtf",wtobj = mod_obj(fba),solverParm=list(CPX_PARAM_EPRHS=1e-07))# tolerance to bounds)))
    wt_flux<-getFluxDist(mtf)
    
    
    model<-constraint_model(model,wt_flux)
    #model<-constraint_model1(model,wt_flux)
    return(model)
}


#============== apply the constraint on the Rich media  ===========================
constraint_richt<-function(model,C_react_pos,C_atoms,rlb,rtype){
 
    fba<-optimizeProb(model ,algorithm="fba1",ReactIndx=C_react_pos, Scoef=C_atoms,rLB=rlb,rtype=rtype) # wt_flux for mtf
    mtf <- optimizeProb(model, algorithm = "mtf_Carbon", ExmetIndx=C_react_pos, cAtoms=C_atoms,rUbLb=rlb, wtobj = mod_obj(fba))
    wt_flux<-getFluxDist(mtf)
    
    
    model<-constraint_model(model,wt_flux)
    
    return(model)
}


#irrevModel<-mod2irrev(Ecoli)
#model<-constraint_minmial(irrevModel)


