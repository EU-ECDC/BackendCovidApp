# Updated version: interpret and check the specifications of each component
# control must contain all arguments, i.e. setControl was used
my.interpretControl <- function (control, stsObj)
{
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  
  Y <- observed(stsObj)
  
  
  ##########################################################################
  ##  get the model specifications for each of the three components
  ##########################################################################
  ar <- control$ar
  ne <- control$ne
  end <- control$end
  
  ## for backwards compatibility with surveillance < 1.8-0, where the ar and ne
  ## components of the control object did not have an offset
  if (is.null(ar$offset)) ar$offset <- 1
  if (is.null(ne$offset)) ne$offset <- 1
  ## for backward compatibility with surveillance < 1.9-0
  if (is.null(ne$normalize)) ne$normalize <- FALSE
  
  # ## create list of offsets of the three components
  Ylagged <- hhh4addon:::weightedSumAR(observed = Y, lag = ar$lag, #BJ
                                       funct_lag = control$funct_lag, par_lag = control$par_lag, min_lag = control$min_lag, max_lag = control$max_lag, #BJ
                                       sum_up = TRUE) #BJ
  if(!is.null(ne$scale)) Ylocal <- local_age(C = ne$scale, Ylagged) else Ylocal <- Ylagged
  
  # This step was just for de-bugging and has been deprecated.
  # Ylagged_detailed <- hhh4addon:::weightedSumAR(observed = Y, lag = ar$lag, #BJ
  #                            funct_lag = control$funct_lag, par_lag = control$par_lag, min_lag = control$min_lag, max_lag = control$max_lag, #BJ
  #                            sum_up = FALSE) #BJ
  
  Ylagged.ne <- hhh4addon:::neOffsetFUN(Y = Y, neweights = ne$weights, scale = ne$scale, normalize = ne$normalize, #BJ
                                           nbmat = neighbourhood(stsObj), data = control$data, lag = control$ne$lag, funct_lag = control$funct_lag, par_lag = control$par_lag, # BJ
                                           min_lag = control$min_lag, max_lag = control$max_lag, sum_up = TRUE, offset = ne$offset)# BJ
  offsets <- list(ar = ar$offset * Ylocal, ne = Ylagged.ne, end = end$offset) #BJ
  
  ## Initial parameter vector 'd' of the neighbourhood weight function
  initial.d <- if (is.list(ne$weights)) ne$weights$initial else numeric(0L)
  dim.d <- length(initial.d)
  names.d <- if (dim.d == 0L) character(0L) else {
    paste0("neweights.", if (is.null(names(initial.d))) {
      if (dim.d==1L) "d" else paste0("d", seq_len(dim.d))
    } else names(initial.d))
  }
  
  ## determine all NA's (FIXME: why do we need this? Why include is.na(Y)?)
  isNA <- is.na(Y)
  if (ar$inModel)
    isNA <- isNA | is.na(offsets[[1L]])
  if (ne$inModel)
    isNA <- isNA | is.na(offsets[[2L]](initial.d))
  
  ## get terms for all components
  all.term <- NULL
  if(ar$isMatrix) stop("matrix-form of 'control$ar$f' is not implemented")
  if(ar$inModel) # ar$f is a formula
    all.term <- cbind(all.term, surveillance:::checkFormula(ar$f, 1, control$data, stsObj))
  if(ne$inModel)
    all.term <- cbind(all.term, surveillance:::checkFormula(ne$f, 2, control$data, stsObj))
  if(end$inModel)
    all.term <- cbind(all.term, surveillance:::checkFormula(end$f,3, control$data, stsObj))
  
  dim.fe <- sum(unlist(all.term["dim.fe",]))
  dim.re.group <- unlist(all.term["dim.re",], use.names=FALSE)
  dim.re <- sum(dim.re.group)
  dim.var <- sum(unlist(all.term["dim.var",]))
  dim.corr <- sum(unlist(all.term["corr",]))
  
  if(dim.corr>0){
    if(dim.var!=dim.corr) stop("Use corr=\'all\' or corr=\'none\' ")
    dim.corr <- switch(dim.corr,0,1,3)
  }
  
  # the vector with dims of the random effects must be equal if they are correlated
  if(length(unique(dim.re.group[dim.re.group>0]))!=1 & dim.corr>0){
    stop("Correlated effects must have same penalty")
  }
  
  n <- c("ar","ne","end")[unlist(all.term["offsetComp",])]
  names.fe <- names.var <- names.re <- character(0L)
  for(i in seq_along(n)){
    .name <- all.term["name",i][[1]]
    names.fe <- c(names.fe, paste(n[i], .name, sep="."))
    if(all.term["random",i][[1]]) {
      names.var <- c(names.var, paste("sd", n[i], .name, sep="."))
      names.re <- c(names.re, paste(n[i], .name, if (.name == "ri(iid)") {
        colnames(stsObj)
      } else {
        seq_len(all.term["dim.re",i][[1]])
      }, sep = "."))
    }
  }
  index.fe <- rep(1:ncol(all.term), times=unlist(all.term["dim.fe",]))
  index.re <- rep(1:ncol(all.term), times=unlist(all.term["dim.re",]))
  
  # poisson or negbin model
  if(identical(control$family, "Poisson")){
    ddistr <- function(y,mu,size){
      dpois(y, lambda=mu, log=TRUE)
    }
    dim.overdisp <- 0L
    index.overdisp <- names.overdisp <- NULL
  } else { # NegBin
    ddistr <- function(y,mu,size){
      dnbinom(y, mu=mu, size=size, log=TRUE)
    }
    ## version that can handle size = Inf (i.e. the Poisson special case):
    ## ddistr <- function (y,mu,size) {
    ##     poisidx <- is.infinite(size)
    ##     res <- y
    ##     res[poisidx] <- dpois(y[poisidx], lambda=mu[poisidx], log=TRUE)
    ##     res[!poisidx] <- dnbinom(y[!poisidx], mu=mu[!poisidx],
    ##                              size=size[!poisidx], log=TRUE)
    ##     res
    ## }
    index.overdisp <- if (is.factor(control$family)) {
      control$family
    } else if (control$family == "NegBinM") {
      factor(colnames(stsObj), levels = colnames(stsObj))
      ## do not sort levels (for consistency with unitSpecific effects)
    } else { # "NegBin1"
      factor(character(nUnits))
    }
    names(index.overdisp) <- colnames(stsObj)
    dim.overdisp <- nlevels(index.overdisp)
    names.overdisp <- if (dim.overdisp == 1L) {
      "-log(overdisp)"
    } else {
      paste0("-log(", paste("overdisp", levels(index.overdisp), sep = "."), ")")
    }
  }
  environment(ddistr) <- getNamespace("stats")  # function is self-contained
  
  # parameter start values from fe() and ri() calls via checkFormula()
  initial <- list(
    fixed = c(unlist(all.term["initial.fe",]),
              initial.d,
              rep.int(2, dim.overdisp)),
    random = as.numeric(unlist(all.term["initial.re",])), # NULL -> numeric(0)
    sd.corr = c(unlist(all.term["initial.var",]),
                rep.int(0, dim.corr))
  )
  # set names of parameter vectors
  names(initial$fixed) <- c(names.fe, names.d, names.overdisp)
  names(initial$random) <- names.re
  names(initial$sd.corr) <- c(names.var, head(paste("corr",1:3,sep="."), dim.corr))
  
  # modify initial values according to the supplied 'start' values
  initial[] <- mapply(
    FUN = function (initial, start, name) {
      if (is.null(start))
        return(initial)
      if (is.null(names(initial)) || is.null(names(start))) {
        if (length(start) == length(initial)) {
          initial[] <- start
        } else {
          stop("initial values in 'control$start$", name,
               "' must be of length ", length(initial))
        }
      } else {
        ## we match by name and silently ignore additional start values
        start <- start[names(start) %in% names(initial)]
        initial[names(start)] <- start
      }
      return(initial)
    },
    initial, control$start[names(initial)], names(initial),
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  
  # Done
  result <- list(response = Y,
                 terms = all.term,
                 nTime = nTime,
                 nUnits = nUnits,
                 nFE = dim.fe,
                 nd = dim.d,
                 nOverdisp = dim.overdisp,
                 nRE = dim.re,
                 rankRE = dim.re.group,
                 nVar = dim.var,
                 nCorr = dim.corr,
                 nSigma = dim.var+dim.corr,
                 nGroups = ncol(all.term),
                 namesFE = names.fe,
                 indexFE = index.fe,
                 indexRE = index.re,
                 initialTheta = c(initial$fixed, initial$random),
                 initialSigma = initial$sd.corr,
                 offset = offsets,
                 family = ddistr,
                 indexPsi = index.overdisp,
                 subset = control$subset,
                 isNA = isNA
  )
  return(result)
}


## Create autoregressive offset: sum contribution of all age groups in region of interest
local_age <- function(C, Y){
  # Extract the number of regions 
  regs <- (sub("[.].*", "", colnames(Y)))
  nb_reg <- length(unique(regs))
  # Compute the within-region matrix of connectivity between age groups
  rownames(C) <- colnames(C) <- colnames(Y)
  for(i in seq_len(nb_reg)){
    # Elements linking different regions are set to 0
    reg_i <- regs[i]
    C[regs == reg_i, regs != reg_i] <- 0 
  }
  # Compute the local offset 
  Ylocal <- Y %*% C
  
  return(Ylocal)
}
