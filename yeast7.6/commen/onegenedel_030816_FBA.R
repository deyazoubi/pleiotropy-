### one gene deletion function with flux reduction 
onegenedel<-function (model,Emodel, geneList,lb = rep(0, length(geneList)),ub = rep(0, length(geneList)), flux_Percent=seq(0,1,0.005) ,checkOptSolObj = FALSE,
fdWT=fd, ...)
{
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (missing(geneList)) {
        if (length(allGenes(model)) < 1) {
            stop("Argument 'geneList' must contain at least one gene!")
        }
        else {
            geneList <- c(1:length(allGenes(model)))
        }
    }
    
    if (is(geneList, "character")) {
        geneList <- match(geneList, allGenes(model))
        if (any(is.na(geneList))) {
            stop("check genelist!")
        }
    }
    if (missing(flux_Percent)) {
        flux_Percent <- 0
    }
    # else
    #  if (is.na(all(match( flux_Percent,seq(0,1,0.10) ))))
    #  {
    #     stop(" check flux Percent they must be one of these values (0, 0.10, 0.20 ,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)  !")
    #}
      
     
   
   
   #fdWT <- .generateFluxdels(model,geneList,1)
   
   #fdWT[["mdlb"]][sapply(fdWT[["mdlb"]], is.null)] <- 0
   #fdWT[["mdub"]][sapply(fdWT[["mdub"]], is.null)] <- 0
   #lbWT1<-fdWT[["mdlb"]]
   #ubWT1<-fdWT[["mdub"]]
   
   # constraint the flux of its corropsponding reactions
   if(flux_Percent==1){
       
        fd <- fdWT
   }else
   {
       fd <- .generateFluxdels(Emodel, geneList,flux_Percent)
   }


# fd <- .generateFluxdels(Emodel, geneList,flux_Percent)
    
    #fd[["mdlb"]][sapply(fd[["mdlb"]], is.null)] <- 0
    #fd[["mdub"]][sapply(fd[["mdub"]], is.null)] <- 0
    #lb1<-fd[["mdlb"]]
    #ub1<-fd[["mdub"]]


#browser()

   #prob<-sysBiolAlg(model,alg='fba1',ReactIndx=fd[["react"]][[1]],Scoef=fdWT[["Scoef"]][[1]],rLB=fd[["RHS"]][[1]],rtype="U")
   #writeProb(problem(prob),fname='test_fba1306.lp')
   
   
   sol <- optimizer1(model = model,react =fd[["react"]] ,ReactIndx=fd[["react"]],Scoef=fdWT[["Scoef"]],rLB=fd[["RHS"]] ,lb = fdWT[["mdlb"]],  ub =  fdWT[["mdub"]] ,rtype=rep("U",length(fd[["RHS"]])),rebuildModel = T, ...)
   
   #sol <- optimizer1(model = model,react =fd[["react"]] ,ReactIndx=fd[["react"]],Scoef=fd[["Scoef"]],rLB=fd[["RHS"]] ,lb = lbWT1,  ub =  ubWT1 ,rtype=rep("U",length(fd[["RHS"]])),rebuildModel = T, ...)

   
  
    optsol <- new("optsol_genedel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
   
      chlb(optsol) <- as.numeric(rep(flux_Percent,length(geneList)))
       chub(optsol) <-as.numeric(rep(flux_Percent,length(geneList)))
   

    
    dels(optsol) <- matrix(allGenes(model)[geneList], ncol = 1)
    fluxdels(optsol) <- fd[["fd"]]
    hasEffect(optsol) <- fd[["heff"]]
    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }
    return(optsol)
}

#optimizeProb(Emodel,ReactIndx=fd[["react"]][[1]],Scoef=fd[["Scoef"]][[1]],rLB=fd[["RHS"]][[1]],alg="fba1")
.generateFluxdels<-function (model, geneList,flux_Percent )
{
    
    message("compute affected fluxes ... ", appendLF = FALSE)
    react <- mapply(sybil::geneDel, geneList, MoreArgs = list(model = model), SIMPLIFY = FALSE)
    heff <- !sapply(react, is.null, simplify = TRUE, USE.NAMES = FALSE)
    fd <- vector(mode = "list", length = length(react))
    mdlb <- vector(mode = "list", length = length(react))
    mdub <- vector(mode = "list", length = length(react))
    
    RHS <- vector(mode = "list", length = length(react))
    Scoef <- vector(mode = "list", length = length(react))
    PScoef <- vector(mode = "list", length = length(react))
    #rev <- vector(mode = "list", length = length(react))
    
    
    fd[heff] <- lapply(react[heff], function(x) react_id(model)[x])
    mdlb[heff] <- lapply(react[heff], function(x) lowbnd (model)[x])
    mdub[heff] <- lapply(react[heff], function(x) uppbnd(model)[x])
    
    PScoef[heff]<-lapply(react[heff], function(x) (lowbnd(model)[x]) + (uppbnd (model)[x]))
    
    Scoef[heff] <- lapply(PScoef[heff], function(x) sign(x))
    
    if(!is.null(mdub[heff]))
    mdub[heff] <- lapply(react[heff], function(x) (uppbnd (model)[x])*flux_Percent)
    mdlb[heff] <- lapply(react[heff], function(x) (lowbnd (model)[x])*flux_Percent)
    
  
    
    RHS[heff]<-lapply(react[heff], function(x) sum((((lowbnd(model)[x])*-flux_Percent) + (uppbnd (model)[x])*flux_Percent)))
    
    
    
    message("OK")
    
    return(list(react = react, heff = heff, fd = fd, mdlb= mdlb, mdub= mdub,RHS=RHS,Scoef=Scoef))
}




