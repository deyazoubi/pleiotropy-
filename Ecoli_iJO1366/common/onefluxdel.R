#### one flux deletion function with flux reduction 
onefluxdel<-function (model, react = c(1:react_num(model)), lb = rep(0, length(react)), 
    ub = rep(0, length(react)), flux_Percent=seq(0,1,0.005), checkOptSolObj = FALSE, ...) 
{
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    creact <- checkReactId(model, react)
    if (!is(creact, "reactId")) {
        stop("check argument react")
    }
    if (missing(flux_Percent)) {
        flux_Percent <- 0
    }
    # else
    #  if (is.na(all(match( flux_Percent,seq(0,1,0.10) ))))
    #  {
    #     stop(" check flux Percent they must be one of these values (0, 0.10, 0.20 ,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)  !")
    #}
    
    creact <- sort(react_pos(creact))
    #browser()
    fd<-.generateFluxdels_ofd(model,creact,flux_Percent)
    lb1<-fd[["mdlb"]]
    ub1<-fd[["mdub"]]
    
    sol <- optimizer(model = model, lb = lb1, ub = ub1, react = as.list(creact), ...)
    optsol <- new("optsol_fluxdel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
    chlb(optsol) <- as.numeric(lb)
    chub(optsol) <- as.numeric(ub)
    dels(optsol) <- matrix(react_id(model)[creact], ncol = 1)
    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }
    return(optsol)
}


.generateFluxdels_ofd<-function (model, ReactList,flux_Percent )
{ 
	
	
	react =as.list(ReactList)
	
    message("compute affected fluxes ... ", appendLF = FALSE)
    #react <- mapply(sybil::geneDel, geneList, MoreArgs = list(model = model), SIMPLIFY = FALSE)
    heff <- !sapply(react, is.null, simplify = TRUE, USE.NAMES = FALSE)
    #fd <- vector(mode = "list", length = length(react))
    mdlb <- vector(mode = "list", length = length(react))
    mdub <- vector(mode = "list", length = length(react))
    
    #RHS <- vector(mode = "list", length = length(react))
    #Scoef <- vector(mode = "list", length = length(react))
    #PScoef <- vector(mode = "list", length = length(react))
    #rev <- vector(mode = "list", length = length(react))
    
    
    #fd[heff] <- lapply(react[heff], function(x) react_id(model)[x])
    mdlb[heff] <- lapply(react[heff], function(x) lowbnd (model)[x])
    mdub[heff] <- lapply(react[heff], function(x) uppbnd(model)[x])
    
    #PScoef[heff]<-lapply(react[heff], function(x) (lowbnd(model)[x]) + (uppbnd (model)[x]))
    
    #Scoef[heff] <- lapply(PScoef[heff], function(x) sign(x))
    
    #if(!is.null(mdub[heff]))
    mdub[heff] <- lapply(react[heff], function(x) (uppbnd (model)[x])*flux_Percent)
    mdlb[heff] <- lapply(react[heff], function(x) (lowbnd (model)[x])*flux_Percent)
    
    # browser()

    #RHS[heff]<-lapply(react[heff], function(x) sum((((lowbnd(model)[x])*-flux_Percent) + (uppbnd (model)[x])*flux_Percent)))
    
    
    #rev[heff] <- lapply(react[heff], function(x) react_rev(model)[x])
     #dd<-lapply(fd[heff], function(x) sapply(strsplit(x,"_"),tail,1))
     #names(Scoef[[1]])<-dd[[1]]
     #which(Scoef[[1]]!=1)
    #sapply(strsplit(fd,"_"),tail,1)
    
    #mdub[sapply(mdub, is.null)] <- 0
    #mdub<-unlist(mdub)

    message("OK")

    return(list(react = react, heff = heff, mdlb= mdlb, mdub= mdub))
}


