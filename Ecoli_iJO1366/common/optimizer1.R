optimizer1<-
  function (model, react, lb = NULL, ub = NULL, obj_coef = NULL, 
            lpdir = NULL, algorithm = SYBIL_SETTINGS("ALGORITHM"), mtfobj = NULL, 
            setToZero = FALSE, rebuildModel = FALSE, fld = "none", prCmd = NA, 
            poCmd = NA, prDIR = NULL, poDIR = NULL, verboseMode = 2,ReactIndx=NULL,Scoef=NULL,rLB=NULL,rtype=NULL, 
            ...) 
  {
   # browser()
    
  
    stopifnot(length(fld) == 1)
    on.exit(expr = {
      if (exists("logObj")) {
        logClose(logObj) <- NA
      }
    })
    logObj <- sybilLog(filename = "", loglevel = -1, verblevel = verboseMode)
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg")
    }
    if (!is(react, "list")) {
      stop("needs an object of class list")
    }
    nObj <- length(react)
    if (!is.null(lb)) {
      stopifnot((length(lb) == nObj || nrow(lb) == nObj), (is(lb, 
                                                              "numeric") || is(lb, "list") || is(lb, "matrix")))
    }
    if (!is.null(ub)) {
      stopifnot((length(ub) == nObj || nrow(ub) == nObj), (is(ub, 
                                                              "numeric") || is(ub, "list") || is(ub, "matrix")))
    }
    if (!is.null(obj_coef)) {
      stopifnot((length(obj_coef) == nObj || nrow(obj_coef) == 
                   nObj), (is(obj_coef, "numeric") || is(obj_coef, "list") || 
                             is(obj_coef, "matrix")))
    }
    if (!is.null(lpdir)) {
      stopifnot(length(lpdir) == nObj, is(lpdir, "character"))
    }
    if (isTRUE(fld)) {
      fdist <- "all"
    }
    else if (identical(fld, FALSE)) {
      fdist <- "none"
    }
    else {
      fdist <- fld
    }
    gEl <- function(el, num, pos) {
      if (is.null(el)) {
        return(el)
      }
      stopifnot(length(pos) == 1)
      if (is.list(el)) {
        ret <- el[[pos]]
      }
      else if (is.matrix(el)) {
        ret <- el[pos, , drop = TRUE]
      }
      else {
        ret <- rep(el[pos], num)
      }
      stopifnot(length(ret) == num)
      return(ret)
    }
    if (algorithm == "mtf") {
      if (fdist == "none") {
        fdist <- "fluxes"
      }
      if (is.null(mtfobj)) {
        lpmod <- sysBiolAlg(model, algorithm = "mtf", react = react, 
                            lb = lb, ub = ub, ...)
      }
      else {
        stopifnot(is(mtfobj, "numeric"), length(mtfobj) == 
                    nObj)
        lpmod <- sysBiolAlg(model, algorithm = "mtf", wtobj = mtfobj, 
                            ...)
      }
    }
    else {
      if (algorithm == "fba1") { 
      #print("sysBiolAlg")
      lpmod <- sysBiolAlg(model, algorithm = 'fba',...)}
      else{
      lpmod <- sysBiolAlg(model, algorithm = algorithm, ...)
      }
    }
    pert <- checkAlgorithm(algorithm, "pert")
    obj <- numeric(nObj)
    mobj <- numeric(nObj)
    ok <- integer(nObj)
    stat <- integer(nObj)
    flux <- switch(fdist, all = {
      Matrix::Matrix(0, nrow = nc(lpmod), ncol = nObj)
    }, fluxes = {
      Matrix::Matrix(0, nrow = length(fldind(lpmod)), ncol = nObj)
    }, {
      NA
    })
    runPrPl <- logical(nObj)
    runPoPl <- logical(nObj)
    runPrPcn <- 1
    runPoPcn <- 1
    if (all(!is.na(prCmd))) {
      do_pr <- TRUE
      prPcmd <- NULL
      runPrP <- .doInRound(prDIR, nObj)
      prPpa <- vector(mode = "list", length = length(runPrP))
      runPrPl[runPrP] <- TRUE
    }
    else {
      do_pr <- FALSE
    }
    if (all(!is.na(poCmd))) {
      do_po <- TRUE
      poPcmd <- NULL
      runPoP <- .doInRound(poDIR, nObj)
      poPpa <- vector(mode = "list", length = length(runPoP))
      runPoPl[runPoP] <- TRUE
    }
    else {
      do_po <- FALSE
    }
    message("calculating ", nObj, " optimizations ... ", appendLF = FALSE)
    if (verboseMode > 1) {
      cat("\n")
    }
    if (verboseMode == 2) {
      progr <- sybil:::.progressBar()
    }
    logOptimizationTH(logObj)
    fi <- fldind(lpmod)
    for (i in 1:nObj) {
      if (verboseMode == 2) {
        progr <- sybil:::.progressBar(i, nObj, progr)
      }
      if (isTRUE(runPrPl[i])) {
        prCmd_tmp <- prCmd
        did_pr <- TRUE
      }
      else {
        prCmd_tmp <- NA
        did_pr <- FALSE
      }
      if (isTRUE(runPoPl[i])) {
        poCmd_tmp <- poCmd
        did_po <- TRUE
      }
      else {
        poCmd_tmp <- NA
        did_po <- FALSE
      }
      
      if (isTRUE(rebuildModel)) {
       #browser()
        sol <- optimizeProb(model, react = react[[i]], lb = gEl(lb, 
                            length(react[[i]]), i), ub = gEl(ub, length(react[[i]]), i), obj_coef = gEl(obj_coef, length(react[[i]]),i), lpdir = lpdir[i],
                            retOptSol = FALSE, prCmd = prCmd_tmp,poCmd = poCmd_tmp, prCil = runPrPcn, poCil = runPoPcn, 
                            algorithm = algorithm,ReactIndx=ReactIndx[[i]],Scoef=Scoef[[i]],rLB=rLB[[i]],rtype=rtype[i], ...)
        
      }
      else {
        if (algorithm == "mtf") {
          changeMaxObj(lpmod, i)
        }
       
        sol <- optimizeProb(lpmod, react = react[[i]], lb = gEl(lb, 
                                                                length(react[[i]]), i), ub = gEl(ub, length(react[[i]]), 
                                                                                                 i), obj_coef = gEl(obj_coef, length(react[[i]]), 
                                                                                                                    i), lpdir = lpdir[i], prCmd = prCmd_tmp, poCmd = poCmd_tmp, 
                            prCil = runPrPcn, poCil = runPoPcn)
      }
      ok[i] <- sol$ok
      stat[i] <- sol$stat
      if (fdist == "none") {
        if (isTRUE(pert)) {
          obj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
        }
        else {
          obj[i] <- sol$obj
        }
      }
      else {
        obj[i] <- sol$obj
        if (fdist == "fluxes") {
          flux[, i] <- sol$fluxes[fi]
        }
        else {
          flux[, i] <- sol$fluxes
        }
      }
      if (isTRUE(did_pr)) {
        if ((runPrPcn == 1) && (is.null(prPcmd))) {
          prPcmd <- cmd(sol$preP)
        }
        prPpa[[runPrPcn]] <- pa(sol$preP)
        runPrPcn <- runPrPcn + 1
        did_pr <- FALSE
      }
      if (isTRUE(did_po)) {
        if ((runPoPcn == 1) && (is.null(poPcmd))) {
          poPcmd <- cmd(sol$postP)
        }
        poPpa[[runPoPcn]] <- pa(sol$postP)
        runPoPcn <- runPoPcn + 1
        did_po <- FALSE
      }
      logOptimization(logObj, sol$ok, sol$stat, obj[i], lpdir[i], 
                      obj_coef[i], react[[i]], i)
      remove(sol)
    }
    message("OK")
    if (fdist == "fluxes") {
      fli <- 1:length(fi)
    }
    else if (fdist == "none") {
      fli <- NA
    }
    else {
      fli <- fi
    }
    if (isTRUE(do_pr)) {
      prAna <- ppProc(prPcmd)
      pa(prAna) <- prPpa
      ind(prAna) <- runPrP
    }
    else {
      prAna <- NULL
    }
    if (isTRUE(do_po)) {
      poAna <- ppProc(poPcmd)
      pa(poAna) <- poPpa
      ind(poAna) <- runPoP
    }
    else {
      poAna <- NULL
    }
    optsol <- list(solver = solver(problem(lpmod)), method = method(problem(lpmod)), 
                   algorithm = algorithm(lpmod), lp_num_cols = nc(lpmod), 
                   lp_num_rows = nr(lpmod), obj = obj, ok = ok, stat = stat, 
                   lp_dir = factor(getObjDir(problem(lpmod))), fldind = fli, 
                   fluxdist = fluxDistribution(flux), prAna = prAna, poAna = poAna, 
                   alg_par = alg_par(lpmod))
    if (isTRUE(setToZero)) {
      do_again <- checkSolStat(stat, solver(problem(lpmod)))
      num_new <- length(do_again)
      optsol[["obj"]][do_again] <- as.numeric(0)
      message("setting ", num_new, " objective values to zero")
      for (i in seq(along = do_again)) {
        logOptimization(logObj, optsol[["ok"]][do_again[i]], 
                        optsol[["stat"]][do_again[i]], 0, lpdir[[do_again[i]]], 
                        obj_coef[[do_again[i]]], react[[do_again[i]]], 
                        do_again[i])
      }
    }
    delProb(problem(lpmod))
    remove(lpmod)
    logFoot(logObj) <- TRUE
    logClose(logObj) <- NA
    return(optsol)
  }

