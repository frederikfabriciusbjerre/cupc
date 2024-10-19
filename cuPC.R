library(pcalg)
library(Rfast)

cu_pc <- function(suffStat, indepTest, alpha, labels, p,
                  fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, 
                  m.max = Inf, u2pd = c("relaxed", "rand", "retry"),
                  skel.method = c("stable", "original", "stable.fast"),
                  conservative = FALSE, maj.rule = FALSE,
                  solve.confl = FALSE, verbose = FALSE) {
  ## Initial Checks
  cl <- match.call()
  if (!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    } else if (p != length(labels)) {
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)") 
    } else {
      message("No need to specify 'p', when 'labels' is given")
    }
  }
  seq_p <- seq_len(p)

  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule) {
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    }

    if (solve.confl) {
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
    }
  }

  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")

  ## Skeleton
  skel <- cu_skeleton(suffStat, indepTest, alpha, labels = labels, NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  skel@call <- cl # so that makes it into result

  ## Orient edges
  if (!conservative && !maj.rule) {
    switch(u2pd,
      "rand" = udag2pdag(skel),
      "retry" = udag2pdagSpecial(skel)$pcObj,
      "relaxed" = udag2pdagRelaxed(skel, verbose = verbose)
    )
  } else { ## u2pd "relaxed" : conservative _or_ maj.rule

    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
      version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose
    )
    udag2pdagRelaxed(pc.$sk,
      verbose = verbose,
      unfVect = pc.$unfTripl
    )
  }
}


cu_skeleton <- function(suffStat, indepTest, alpha, labels, p, m.max = Inf, NAdelete = TRUE, verbose = FALSE) {
  cl <- match.call()
  if (!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    } else if (p != length(labels)) {
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    } else {
      message("No need to specify 'p', when 'labels' is given")
    }
  }

  seq_p <- seq_len(p)
  # pval <- NULL
  # Convert SepsetMatrix to sepset
  sepset <- lapply(seq_p, function(.) vector("list", p)) # a list of lists [p x p]
  # save maximal p value
  pMax <- matrix(0, nrow = p, ncol = p)
  number_of_levels <- 50
  suffStat_n <- suffStat[length(suffStat)][[1]]
  threshold <- matrix(0, nrow = 1, ncol = number_of_levels)
  for (i in 0:(min(number_of_levels, suffStat_n - 3) - 1)) {
    threshold[i] <- abs(qnorm((alpha / 2), mean = 0, sd = 1) / sqrt(suffStat_n - i - 3))
  }
  Z_list <- vector("list", M)
  sepset_list <- vector("list", M)  
  ord_list <- numeric(M)            
  for (m in 1:M) {
    C_m <- suffStat[[m]] # Covariance matrix for the m-th imputation
    # Initialize the adjacency matrix G for this imputation
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
    G <- G * 1
    ord <- 0
    max_level <- if (is.infinite(m.max)) 14 else m.max
    
    # Initialize pMax and sepsetMatrix
    pMax <- matrix(0, nrow = p, ncol = p)
    sepsetMatrix <- matrix(-1, nrow = p * p, ncol = 14)
    
    # Load the CUDA library if not already loaded
    if (!is.loaded("Skeleton")) {
      dyn.load("Skeleton.so")
    }
    
    # Call the existing CUDA function to compute Z-statistics for this imputation
    z <- .C("Skeleton",
            C = as.double(C_m),
            p = as.integer(p),
            G = as.integer(G),
            Th = as.double(threshold),
            l = as.integer(ord),
            max_level = as.integer(max_level),
            pmax = as.double(pMax),
            sepsetmat = as.integer(sepsetMatrix))
    
    # Extract ord, G, pMax, sepsetMatrix
    ord_m <- z$l
    ord_list[m] <- ord_m
    G_m <- matrix(z$G, nrow = p, ncol = p)
    
    # Extract Z-statistics from pMax
    Z_m <- matrix(z$pmax, nrow = p, ncol = p)
    Z_m[Z_m == -100000] <- -Inf  # Replace placeholder values with -Inf
    
    # Store the Z-statistics for this imputation
    Z_list[[m]] <- Z_m
    
    # Process sepsetMatrix to create sepset for this imputation
    sepset_m <- lapply(seq_p, function(.) vector("list", p))
    if (ord_m <= 14) {
      sepsetMatrix_m <- t(matrix(z$sepsetmat, nrow = 14, ncol = p^2))
      index_of_cut_edge <- which(sepsetMatrix_m != -1, arr.ind = TRUE)
      for (idx in seq_len(nrow(index_of_cut_edge))) {
        i <- index_of_cut_edge[idx, 1]
        j <- index_of_cut_edge[idx, 2]
        x <- ((i - 1) %% p) + 1
        y <- ((i - 1) %/% p) + 1
        # Extract the separating set
        sep_set <- sepsetMatrix_m[i, 1:ord_m]
        sepset_m[[x]][[y]] <- sepset_m[[y]][[x]] <- sep_set[sep_set != -1]
      }
    } else {
      # TODO: Handle cases where ord > 14
    }
    sepset_list[[m]] <- sepset_m
  }
  print(ord)
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G, "graphNEL")
    }
  ## final object
  new("pcAlgo",
    graph = Gobject, call = cl, n = integer(0),
    max.ord = as.integer(ord - 1), n.edgetests = 0,
    sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1)
  )
} ## end{ skeleton }

# copied and changed from micd github
getSuff <- function(X) {
  if (!(mice::is.mids(X) | is.list(X))) {
    stop("data is neither a list nor a mids object")
  }
  if (inherits(X, "mids")) {
    X <- mice::complete(X, action="all")
    if(length(which.is(X, "factor") > 0)) {
      stop("data must be all numeric.")
    }
  }
  C <- lapply(X, stats::cor)
  n <- nrow(X[[1]])
  c(C, n)
}