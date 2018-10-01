doublegenedel<-function (model,Emodel, geneList1, geneList2, lb = NULL, ub = NULL, 
    allComb = FALSE, exLethal = TRUE, tol = SYBIL_SETTINGS("TOLERANCE"), 
    checkOptSolObj = FALSE, flux_Percent1=seq(0,1,0.01), flux_Percent2=seq(0,1,0.01),fdWT1=fd1,fdWT2=fd2, ...)
{
   
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (missing(geneList1)) {
        geneList1 <- allGenes(model)
    }
    else {
        if (is.na(all(match(geneList1, allGenes(model))))) {
            stop("check genelist1!")
        }
    }
    if (missing(geneList2)) {
        geneList2 <- allGenes(model)
    }
    else {
        if (is.na(all(match(geneList2, allGenes(model))))) {
            stop("check genelist2!")
        }
    }
    if (length(geneList1) < 1) {
        stop("Argument 'geneList1' must contain at least one gene!")
    }
    if (length(geneList2) < 1) {
        stop("Argument 'geneList2' must contain at least one gene!")
    }
    
    if (missing(flux_Percent1)) {
        flux_Percent1 <- 0
    }
    if (missing(flux_Percent2)) {
        flux_Percent2 <- 0
    }
   
    
    geneList1 <- match(geneList1, allGenes(model))
    geneList2 <- match(geneList2, allGenes(model))
    num_geneList1 <- length(geneList1)
    num_geneList2 <- length(geneList2)
    if ((!isTRUE(allComb)) && (num_geneList1 != num_geneList2)) {
        stop("geneList1 and geneList2 must have the same length")
    }
    num_genes <- length(allGenes(model))
    if (isTRUE(exLethal)) {
        unGenes <- sort(unique(c(geneList1, geneList2)))
        ca <- match.call()
        if ("solver" %in% names(ca)) {
            testslv <- tryCatch(eval(parse(text = ca["solver"])), 
                error = function(e) e)
            if (is(testslv, "simpleError")) {
                slv <- as.character(ca["solver"])
            }
            else {
                slv <- as.character(testslv)
            }
        }
        else {
            slv <- SYBIL_SETTINGS("SOLVER")
        }
        wtobj <- optimizeProb(model, retOptSol = FALSE, algorithm = "fba", 
            solver = slv, lpdir = "max")
        unSol <- oneGeneDel(model, geneList = unGenes, solver = slv, 
            lpdir = rep("max", length(unGenes)), fld = "none", 
            algorithm = "fba")
        letid <- lethal(unSol, wtobj$obj)
        statNok <- checkStat(unSol[letid])
        lethal <- letid
        lethalGeneIds <- unGenes[lethal]
        l1 <- which(geneList1 %in% lethalGeneIds)
        l2 <- which(geneList2 %in% lethalGeneIds)
        if (isTRUE(allComb)) {
            if (length(l1) > 0) {
                geneList1 <- geneList1[-l1]
            }
            if (length(l2) > 0) {
                geneList2 <- geneList2[-l2]
            }
        }
        else {
            l12 <- sort(unique(c(l1, l2)))
            if (length(l12) > 0) {
                geneList1 <- geneList1[-l12]
                geneList2 <- geneList2[-l12]
            }
        }
        remove(unSol)
    }
    else {
        lethalGeneIds <- as.integer(NA)
    }
    if (isTRUE(allComb)) {
        tmpMAT <- upper.tri(matrix(nrow = num_genes, ncol = num_genes), 
            diag = FALSE)
        tmpDIFF <- setdiff(geneList2, geneList1)
        if (length(tmpDIFF) != 0) {
            tmpMAT[, tmpDIFF] <- TRUE
        }
        tmpDIFF <- setdiff(geneList1, geneList2)
        if (length(tmpDIFF) != 0) {
            tmpMAT[tmpDIFF, ] <- TRUE
        }
        geneList1 <- unique(geneList1)
        geneList2 <- unique(geneList2)
        tmpMAT <- tmpMAT[geneList1, geneList2, drop = FALSE]
        num_opt <- sum(tmpMAT == TRUE)
        rownames(tmpMAT) <- geneList1
        colnames(tmpMAT) <- geneList2
        deletions <- which(tmpMAT == TRUE, arr.ind = TRUE)
        kogenesID <- cbind(geneList1[deletions[, "row"]], geneList2[deletions[, 
            "col"]])
    }
    else {
        kogenesID <- cbind(geneList1, geneList2)
    }
    
  
     
    
    kogenes <- lapply(seq_len(nrow(kogenesID)), function(x) kogenesID[x, 
        ])
        #fd <- .generateFluxdels(model, kogenes)
        
        # constraint the flux of its corropsponding reactions
        
     
        
        
        if(flux_Percent1==1){
            
            fd1 <- fdWT1
            
        }else if(is.null(fdWT1[["react"]][[1]])){
		
			fd1 <- fdWT1
        }else if(length(intersect(fdWT1[["react"]][[1]],fdWT2[["react"]][[1]])>0)){		
		
		fd1 <- .generateFluxdels(Emodel, kogenes[[1]][1],min(flux_Percent1,flux_Percent2))
        
        }else{
            fd1 <- .generateFluxdels(Emodel, kogenes[[1]][1],flux_Percent1)
        }
        
        
        if(flux_Percent2==1){
            
            fd2 <- fdWT2
        }else if(is.null(fdWT2[["react"]][[1]])){
		
			fd2 <- fdWT2
        }else if(length(intersect(fdWT1[["react"]][[1]],fdWT2[["react"]][[1]])>0)){		
		
		fd2 <- .generateFluxdels(Emodel, kogenes[[1]][2],min(flux_Percent1,flux_Percent2))
        
        }else{
            fd2 <- .generateFluxdels(Emodel, kogenes[[1]][2],flux_Percent2)
        }
  
    
    
    
    
    ###return to this point
   
       if (is.null(fdWT1[["mdlb"]])) {
        lb <- rep(0, length(kogenes))
    }
    else {
        if (length(fdWT1[["mdlb"]]) != length(kogenes)) {
            stop("lb must be of length ", length(kogenes))
        }
    }
    if (is.null(fdWT1[["mdub"]])) {
        ub <- rep(0, length(kogenes))
    }
    else {
        if (length(fdWT1[["mdub"]]) != length(kogenes)) {
            stop("ub must be of length ", length(kogenes))
        }
    }
    
    
   
    # first case NULL NULL  Gs		
    if(is.null(fd1[["react"]][[1]]) & is.null(fd2[["react"]][[1]])){
		
		           
     fdWT <- .generateFluxdels(model,kogenes,1)#rep(1,length((fd[["react"]][[1]])))
        
     fd <- .generateFluxdels(Emodel,kogenes,min(sum(flux_Percent1,flux_Percent2),1))# Min(sum(Flux%1, Flux%2),1)
          
	  fd[["mdlb"]][sapply(fd[["mdlb"]], is.null)] <- 0
	  fd[["mdub"]][sapply(fd[["mdub"]], is.null)] <- 0
		 
	  if(is.null(fd[["react"]][[1]])){
		  
		  prob<-sysBiolAlg(model)
		  opt<-optimizeProb(prob)
		  return( list(opt,F))
		 
		 }else{
			 
		 ec <- list(
			react= fd[["react"]],
			x=list(fdWT[["Scoef"]][[1]]), ub =  fd[["RHS"]][[1]],
			rtype="U")
			
		sol <- optimizer(model = model,react = fd[["react"]],easyConstraint=ec, ...)
			#browser()
			flag<-T
		
	}
	
	# case two NULL Gs 	
	}else if(is.null(fd1[["react"]][[1]])){
		#G1 <-NULL
		 fdWT <- .generateFluxdels(model,kogenes,1)
		 
		if(length(fdWT[["react"]][[1]])> length(fd2[["react"]][[1]]) ){
			# G1 NULL G2 genes @ G1G2 genes+1
			fd1 <- .generateFluxdels(Emodel, kogenes[[1]][1],min(sum(flux_Percent1,flux_Percent2),1))
			fd2 <- .generateFluxdels(Emodel, kogenes[[1]][2],min(sum(flux_Percent1,flux_Percent2),1))
			
			ec <- list(
			react=list( fd1[["react"]][[1]], fd2[["react"]][[1]]),
			x=list(fdWT1[["Scoef"]][[1]],fdWT2[["Scoef"]][[1]] ), ub =  c(fd1[["RHS"]][[1]],fd2[["RHS"]][[1]]),
			rtype=c("U","U"))
    
    prob<-sysBiolAlg(model,easyConstraint=ec, ...)
    opt<-optimizeProb(prob)
     flag<-F
			
			
		}else{
		ec <- list(
			react= fd2[["react"]],
			x=list(fdWT2[["Scoef"]][[1]]), ub =  fd2[["RHS"]][[1]],
			rtype="U")
			
		sol <- optimizer(model = model,react = fd2[["react"]],easyConstraint=ec, ...)
		flag<-T
	}
		
	}else if(is.null(fd2[["react"]][[1]])){
	# G2 NULL G1 genes @ G1G2 genes+1
	fdWT <- .generateFluxdels(model,kogenes,1)
	if(length(fdWT[["react"]][[1]]) > length(fd1[["react"]][[1]]) ){
			
			fd1 <- .generateFluxdels(Emodel, kogenes[[1]][1],min(sum(flux_Percent1,flux_Percent2),1))
			fd2 <- .generateFluxdels(Emodel, kogenes[[1]][2],min(sum(flux_Percent1,flux_Percent2),1))
			
			ec <- list(
			react=list( fd1[["react"]][[1]], fd2[["react"]][[1]]),
			x=list(fdWT1[["Scoef"]][[1]],fdWT2[["Scoef"]][[1]] ), ub =  c(fd1[["RHS"]][[1]],fd2[["RHS"]][[1]]),
			rtype=c("U","U"))
    
    prob<-sysBiolAlg(model,easyConstraint=ec, ...)
    opt<-optimizeProb(prob)
     flag<-F
			
			
		}else{
			
			
	ec <- list(
			react= fd1[["react"]],
			x=list(fdWT1[["Scoef"]][[1]]), ub =  fd1[["RHS"]][[1]],
			rtype="U")
			
		sol <- optimizer(model = model,react = fd1[["react"]],easyConstraint=ec, ...)
		flag<-T
		}
		
	}else if(length(intersect(fd1[["react"]][[1]],fd2[["react"]][[1]])>0)){
		
		ec <- list(
			react=list( fd1[["react"]][[1]], fd2[["react"]][[1]]),
			x=list(fdWT1[["Scoef"]][[1]],fdWT2[["Scoef"]][[1]] ), ub =  c(fd1[["RHS"]][[1]],fd2[["RHS"]][[1]]),
			rtype=c("U","U"))
    
    prob<-sysBiolAlg(model,easyConstraint=ec, ...)
    opt<-optimizeProb(prob)
     flag<-F

	}
	else{
	 #browser()
    ec <- list(
			react=list( fd1[["react"]][[1]], fd2[["react"]][[1]]),
			x=list(fdWT1[["Scoef"]][[1]],fdWT2[["Scoef"]][[1]] ), ub =  c(fd1[["RHS"]][[1]],fd2[["RHS"]][[1]]),
			rtype=c("U","U"))
    
    prob<-sysBiolAlg(model,easyConstraint=ec, ...)
    opt<-optimizeProb(prob)
     flag<-F
     
 }
       
       
       if(flag){
       optsol <- new("optsol_genedel")
		opt <- makeOptsolMO(model, sol)
		as(optsol, "optsol_optimizeProb") <- opt
		chlb(optsol) <- as.numeric(lb)
		chub(optsol) <- as.numeric(ub)
		dels(optsol) <- matrix(allGenes(model)[kogenesID], ncol = 2)
		#fluxdels(optsol) <- fd[["fd"]]
		#hasEffect(optsol) <- fd[["heff"]]
		
		
			if (isTRUE(checkOptSolObj)) {
				checkOptSol(optsol, onlywarn = TRUE)
			}
    
    
    
    return(list(optsol,T))
    
    
}else{
	
	return( list(opt,F))
	
}
    
}

   
   

  
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


