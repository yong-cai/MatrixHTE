# Functions for various bounds on treatment effect on matrices.
require(ggplot2)
# ==============================================================================
# Content
# ------------------------------------------------------------------------------
#
# ==============================================================================

# Smoothed indicator functions
# Currently only implemented with normal CDF
# If bandwidth not provided, use n^{-1/5}
ind_smooth <- function(mat, h=NULL){
  # mat <- Treat.Mat.Sym
  if(is.null(h)){
    h <- dim(mat)[1]^(-1/5)
  }
  return(1-pnorm(mat/h))
}

# Take two matrices of different sizes and return them blown up to equal sizes
# Accepts both symmetric and bipartite networks
joint_blowup <- function(Treat.Mat, Control.Mat){

  treat.n <- dim(Treat.Mat)
  control.n <- dim(Control.Mat)

  # Compute scaling factors and scale the matrices
  treat.scale.1 <- pracma::Lcm(treat.n[1], control.n[1])/treat.n[1]
  control.scale.1 <- pracma::Lcm(treat.n[1], control.n[1])/control.n[1]
  treat.scale.2 <- pracma::Lcm(treat.n[2], control.n[2])/treat.n[2]
  control.scale.2 <- pracma::Lcm(treat.n[2], control.n[2])/control.n[2]

  treat.out <- kronecker(Treat.Mat, matrix(1,treat.scale.1,treat.scale.2))
  control.out <- kronecker(Control.Mat, matrix(1,control.scale.1,control.scale.2))

  return(list(treat.out, control.out))

}

# Symmetrize a non-square matrix square
matrix_symmetrize <- function(A, force = FALSE){
  m <- dim(A)[1]
  n <- dim(A)[2]

  if(m == n & force == FALSE){
    sym <- A
  } else{
    sym <- rbind(cbind(matrix(0, m, m), A),
                 cbind(t(A),matrix(0, n, n)))
  }
  return(sym)
}

joint_blowup_eigen <- function(eig1, eig0, method = "discrete"){

  n.treat <- length(eig1)
  n.control <- length(eig0)

  if(method == "discrete"){
    if(n.treat >= n.control){
      eig1.blowup <- n.treat/n.control*sort(c(eig0, rep(0, n.treat-n.control)), decreasing = TRUE)
      eig0.blowup <- n.control/n.treat*sort(eig1[order(abs(eig1),decreasing = TRUE)][1:n.control], decreasing = TRUE)
    } else {
      eig1.blowup <- n.treat/n.control*sort(eig0[order(abs(eig0),decreasing = TRUE)][1:n.treat], decreasing = TRUE)
      eig0.blowup <- n.control/n.treat*sort(c(eig1, rep(0, n.control-n.treat)), decreasing = TRUE)
    }
  } else if(method == "continuous"){
    eig1.blowup <- (n.treat/n.control)*quantile(eig0, probs = seq(1, 0, length.out = n.treat))
    eig0.blowup <- (n.control/n.treat)*quantile(eig1, probs = seq(1, 0, length.out = n.control))
  } else {
    stop("Method must be 'discrete' or 'continuous'.")
  }

  return(list(eig1.blowup, eig0.blowup))
}

# impute fixed effects.
# for each vector, replace every entry with the corresponding quantile in the other vector
fe.cf <- function(S.means, A.means){
  S.cf <- as.numeric(quantile(A.means, probs = rank(S.means)/length(S.means)))
  A.cf <- as.numeric(quantile(S.means, probs = rank(A.means)/length(A.means)))
  return(list(S.cf, A.cf))
}

# Check if a matrix is square
isSquare <- function(mat){
  dimensions <- dim(mat)
  if(length(dimensions) != 2){
    stop("Input is not of dimension 2")
  }
  return(dimensions[1] == dimensions[2])
}

# Decompose Matrix into Row FE, Column FE and Demeaned Component.
# Assign half the mean to rows and half to columns
matrix_demean <- function(mat){

  mat.dim <- dim(mat)
  m <- mat.dim[1]
  n <- mat.dim[2]
  row.means <- apply(mat, 1, mean) - mean(mat)/2
  col.means <- apply(mat, 2, mean) - mean(mat)/2
  demeaned.mat <- mat - matrix(row.means, m, n, byrow = FALSE) -
    matrix(col.means, m, n, byrow = TRUE) #+ mean(mat)
  return(list("demeaned.mat" = demeaned.mat, "row.means" = row.means, "col.means" = col.means))
}


# Vector-Frechet Bounds for Symmetric and Asymmetric Matrices
# Implemented in terms of eigen/singular values
vector_frechet <- function(Treat.Mat, Control.Mat){

  blowup <- joint_blowup(Treat.Mat, Control.Mat)
  Treat.Mat <- blowup[[1]]
  Control.Mat <- blowup[[2]]

  Treat.ev <- svd(Treat.Mat)$d
  Control.ev <- svd(Control.Mat)$d

  n <- dim(Treat.Mat)[1] # only for symmetric matrices
  m <- dim(Treat.Mat)[2]

  VF_lower <- max(sum(Treat.ev^2) + sum(Control.ev^2) - n*m, 0)/(n*m)
  VF_upper <- min(sum(Treat.ev^2), sum(Control.ev^2))/(n*m)

  return(c(VF_lower, VF_upper))
}

# Eigenvalue Bounds for Symmetric and Asymmetric Matrices
eigenvalue_bounds <- function(Treat.Mat, Control.Mat){

  # Treat.Mat <- S.entry; Control.Mat <- A.entry

  blowup <- joint_blowup(Treat.Mat, Control.Mat)
  Treat.Mat <- blowup[[1]]
  Control.Mat <- blowup[[2]]

  denom <- dim(Treat.Mat)
  denom_adjust <- 1 # factor for adjusting value of denominator. 1 if matrices are symmetric.

  if(denom[1] != denom[2] || !isSymmetric(Treat.Mat) || !isSymmetric(Control.Mat)){
    Treat.Mat <- matrix_symmetrize(Treat.Mat, force = TRUE)
    Control.Mat <- matrix_symmetrize(Control.Mat, force = TRUE)
    denom_adjust <- 2 # 2 since symmetrization doubles the number of non-zero entries
  }

  # Matrices are symmetric at this point. Fix numerical issues causing asymmetry.
  Treat.Mat <- (Treat.Mat + t(Treat.Mat))/2
  Control.Mat <- (Control.Mat + t(Control.Mat))/2

  eigenM = eigen(Treat.Mat, only.values = TRUE)$values
  eigenF = eigen(Control.Mat, only.values = TRUE)$values
  #these guys are already sorted

  eF_lower = (eigenM%*%rev(eigenF))/(denom_adjust*denom[1]*denom[2])
  eF_upper = (eigenM%*%eigenF)/(denom_adjust*denom[1]*denom[2])

  return(c(eF_lower,eF_upper))
}



# Matrix Frechet Hoeffding Bounds
# With option to correct for row and column heterogeneity
#' Bounds for Distribution of Potential Outcomes
#'
#' This function implements the bounds for the distribution of potential outcomes as described in Proposition 1 of Auerbach and Cai (2022).
#' @param Treat.Mat Matrix of outcomes for treated units
#' @param Control.Mat Matrix of outcomes for control units
#' @param y1 y1 in F(y1, y0).
#' @param y0 y0 in F(y1, y0).
#' @param smooth Logical indicating whether the smoothing should be performed. Smoothing is recommended if potential outcomes are estimated.
#' @param bandwidth Bandwidth to be used for smoothing
#' @param hc Logical indicating whether or not to incorporate row and column heterogeneity. See Section 5.2 of Auerbach and Cai (2022).
#' @return Lower and upper bounds for the DPO.
#' @export
#' @examples
#' # Generate data
#' Treat.Mat <- matrix(runif(200) > 0.8, 10, 20)*1
#' Control.Mat <- matrix(runif(100) > 0.8, 10, 10)*1
#'
#' # P(Treat.Mat <= 0.5, Control.Mat <= 0.5)
#' dpo_bounds(Treat.Mat, Control.Mat, 0.5, 0.5)
#'
#' # P(Treat.Mat == 1, Control.Mat == 1)
#' dpo_bounds(-Treat.Mat, -Control.Mat, -0.5, -0.5) # input matrices are binary
dpo_bounds <- function(Treat.Mat, Control.Mat, y1, y0, smooth = FALSE, bandwidth = NULL, hc = FALSE){
  blowup <- joint_blowup(Treat.Mat, Control.Mat)
  Treat.Mat <- blowup[[1]]
  Control.Mat <- blowup[[2]]

  denom <- dim(Treat.Mat)

  if(smooth & is.null(bandwidth)){
    bandwidth <- (denom[1]*denom[2])^(-1/5)
  }

  if(smooth){
    Treat.Mat <- ind_smooth(Treat.Mat - y1, bandwidth)
    Control.Mat <- ind_smooth(Control.Mat - y0, bandwidth)
  } else {
    Treat.Mat <- (Treat.Mat <= y1)
    Control.Mat <- (Control.Mat <= y0)
  }

  if(hc){
    Treat.demean <- matrix_demean(Treat.Mat)
    Control.demean <- matrix_demean(Control.Mat)
    Treat.Mat <- Treat.demean$demeaned.mat
    Control.Mat <- Control.demean$demeaned.mat
  }

  result <- matrix(NA, 2, 2)
  result[1, ] <- vector_frechet(Treat.Mat, Control.Mat)
  result[2, ] <- eigenvalue_bounds(Treat.Mat, Control.Mat)

  if(hc){ # Add back the fixed effects
    Treat.row.means <- sort(Treat.demean$row.means, decreasing = FALSE)
    Treat.col.means <- sort(Treat.demean$col.means, decreasing = FALSE)
    Control.row.means <- sort(Control.demean$row.means, decreasing = FALSE)
    Control.col.means <- sort(Control.demean$col.means, decreasing = FALSE)

    lower <- mean(Treat.row.means*rev(Control.row.means)) + mean(Treat.col.means*rev(Control.col.means)) +
      mean(Treat.row.means)*mean(Control.col.means) + mean(Treat.col.means)*mean(Control.row.means)
    upper <- mean(Treat.row.means*(Control.row.means)) + mean(Treat.col.means*(Control.col.means)) +
      mean(Treat.row.means)*mean(Control.col.means) + mean(Treat.col.means)*mean(Control.row.means)

    result <- result + cbind(rep(lower, 2), rep(upper, 2))
  }

  result <- rbind(result, c(0, 1))
  colnames(result) <- c("lower", "upper")
  rownames(result) <- c("vf", "evb", "trivial")

  return(list("bounds" = c(max(result[, 1]), min(result[,2]))))
}

# Makarov Bounds for Symmetric and Bipartite Matrices
#' Bounds for Distribution of Treatment Effects at a Point
#'
#' This function computes the bounds for the distribution of treatment effects (DTE) at a given point. Specifically, it computes bounds for \eqn{\Delta(y)} for a specified \eqn{y}
#' @param Treat.Mat Matrix of outcomes for treated units
#' @param Control.Mat Matrix of outcomes for control units
#' @param y Point at which to bound DTE
#' @param smooth Logical indicating whether the smoothing should be performed. Smoothing is recommended if potential outcomes are estimated.
#' @param bandwidth Bandwidth to be used for smoothing
#' @param n.grid Number of grid points to search over for evaluating the bounds at each point.
#' @param hc Logical indicating whether or not to incorporate row and column heterogeneity. See Section 5.2 of Auerbach and Cai (2022).
#' @return A list of bounds and the components.\tabular{ll}{
#'    \code{bounds} \tab Lower and upper Makarov bounds for the DTE at \eqn{y} \cr
#'    \tab \cr
#'    \code{breakdown} \tab Breakdown of the bounds into the vectorized Frechet-Hoeffding (vf) and eigenvalue-based (evb) components. Trivial bounds for probability are 0 and 1. \cr
#' }
#' @keywords internal
#' @noRd
#' @examples
#' # Generate data
#' Treat.Mat <- matrix(runif(200) > 0.8, 10, 20)*1
#' Control.Mat <- matrix(runif(100) > 0.8, 10, 10)*1
#'
#' # Generate DTE and plot it
#' bounds <- matrix_makarov(Treat.Mat, Control.Mat, 0.5)
matrix_makarov <- function(Treat.Mat, Control.Mat, y, smooth = FALSE, bandwidth = NULL, n.grid = 20, hc = FALSE){

  blowup <- joint_blowup(Treat.Mat, Control.Mat)
  Treat.Mat <- blowup[[1]]
  Control.Mat <- blowup[[2]]

  denom <- dim(Treat.Mat)

  if(smooth & is.null(bandwidth)){
    bandwidth <- (denom[1]*denom[2])^(-1/5)
  }
  special.case <- c(min(Treat.Mat) - 1e-5,  # Set Treat.i = 0
                    max(Treat.Mat) + 1e-5, # Set Treat.i = 1
                    y+max(Control.Mat) + 1e-5, # Set Control.i = 1
                    y+min(Control.Mat) - 1e-5)# Set Control.i = 0

  action.region <- c(Treat.Mat, y + Control.Mat)

  if(smooth){
    special.case <- special.case + 5*c(-bandwidth, bandwidth, bandwidth, -bandwidth)
  }

  action.region <- action.region[action.region >= special.case[1] & action.region <= special.case[3]]
  if(length(action.region) > 0){
    grid <- c(special.case, quantile(action.region, probs = seq(0, 1, length.out = n.grid)))
  } else {
    grid <- special.case
  }

  lower.bound <- matrix(NA, 2, length(grid))
  upper.bound <- matrix(NA, 2, length(grid))

  for(i in 1:length(grid)){

    y1 <- grid[i]
    y0 <- y-y1

    if(smooth){
      Treat.i <- ind_smooth(Treat.Mat - y1, bandwidth)
      Control.i <- ind_smooth(Control.Mat + y0, bandwidth)
    } else {
      Treat.i <- (Treat.Mat <= y1)
      Control.i <- (-Control.Mat >= y0)
    }

    Treat.ev <- svd(Treat.i)$d/sqrt(denom[1]*denom[2])
    Control.ev <- svd(Control.i)$d/sqrt(denom[1]*denom[2])

    if(hc){
      Treat.demean <- matrix_demean(Treat.i)
      Control.demean <- matrix_demean(Control.i)
      Treat.i <- Treat.demean$demeaned.mat
      Control.i <- Control.demean$demeaned.mat
    }

    vf.bounds <- vector_frechet(Treat.i, Control.i)   # accepts both symmetric and asymmetric matrices
    evb.bounds <- eigenvalue_bounds(Treat.i, Control.i)   # accepts both symmetric and asymmetric matrices

    if(hc){
      Treat.row.means <- sort(Treat.demean$row.means, decreasing = FALSE)
      Treat.col.means <- sort(Treat.demean$col.means, decreasing = FALSE)
      Control.row.means <- sort(Control.demean$row.means, decreasing = FALSE)
      Control.col.means <- sort(Control.demean$col.means, decreasing = FALSE)

      lower <- mean(Treat.row.means*rev(Control.row.means)) + mean(Treat.col.means*rev(Control.col.means)) +
        mean(Treat.row.means)*mean(Control.col.means) + mean(Treat.col.means)*mean(Control.row.means)
      upper <- mean(Treat.row.means*(Control.row.means)) + mean(Treat.col.means*(Control.col.means)) +
        mean(Treat.row.means)*mean(Control.col.means) + mean(Treat.col.means)*mean(Control.row.means)

      vf.bounds <- vf.bounds + c(lower, upper)
      evb.bounds <- evb.bounds + c(lower, upper)

    }

    lower.bound[1, i] <- sum(Treat.ev^2) - vf.bounds[2]
    lower.bound[2, i] <- sum(Treat.ev^2) - evb.bounds[2]

    upper.bound[1, i] <- 1 - sum(Control.ev^2) + vf.bounds[2]
    upper.bound[2, i] <- 1 - sum(Control.ev^2) + evb.bounds[2]

  }

  lower.best <- apply(lower.bound, 1, max)
  upper.best <- apply(upper.bound, 1, min)

  bound_types <- cbind(c(lower.best, 0), c(upper.best, 1))
  rownames(bound_types) <- c("vf", "evb", "trivial")
  colnames(bound_types) <- c("lower", "upper")

  out <- c(max(bound_types[, 1]), min(bound_types[, 2]))

  return(list("bounds" = out, "breakdown" = bound_types))

}


# DTE for Symmetric and Bipartite Matrices
# Wrapper for matrix_makarov
#' Bounds for Distribution of Treatment Effects
#'
#' This function computes the bounds for the distribution of treatment effects (DTE) along a sequence of points. Specifically, it computes bounds for \eqn{\Delta(y)} for \eqn{y \in [start, end]}. The output of this function can be plotted directly using \code{plot()}.
#' @param Treat.Mat Matrix of outcomes for treated units
#' @param Control.Mat Matrix of outcomes for control units
#' @param start Starting value for evaluating the DTE
#' @param end End value for evaluating the DTE
#' @param dte.n.grid Number of points between start and end at which the CDF is bounded.
#' @param smooth Logical indicating whether the smoothing should be performed. Smoothing is recommended if potential outcomes are estimated.
#' @param bandwidth Bandwidth to be used for smoothing
#' @param n.grid Number of grid points to search over for evaluating the bounds at each point.
#' @param hc Logical indicating whether or not to incorporate row and column heterogeneity. See Section 5.2 of Auerbach and Cai (2022).
#' @param eval.points Vector indicating additional points at which DTE should be explicitly evaluated.
#' @return Lower and upper bounds for the DTE over the grid points.
#' @export
#' @examples
#' # Generate data
#' Treat.Mat <- matrix(runif(200) > 0.8, 10, 20)*1
#' Control.Mat <- matrix(runif(100) > 0.8, 10, 10)*1
#'
#' # Generate DTE and plot it
#' bounds <- dte_bounds(Treat.Mat, Control.Mat, -5, 5, smooth = TRUE)
#' plot(bounds)
dte_bounds <- function(Treat.Mat, Control.Mat, start, end, dte.n.grid = 30,
                smooth = FALSE, bandwidth = NULL, n.grid = 20, hc = FALSE, eval.points = NULL){
  # start <- -30; end <- 15; dte.n.grid <- 30; smooth <- TRUE; n.grid <- 20; hc = FALSE
  dte.grid <- seq(start, end, length.out = dte.n.grid)

  if(!is.null(eval.points)){
    dte.grid <- c(dte.grid, eval.points)
    dte.n.grid <- length(dte.grid)
  }

  out <- matrix(NA, dte.n.grid, 3)


  for(i in 1:dte.n.grid){
    y <- dte.grid[i]
    bounds <- matrix_makarov(Treat.Mat, Control.Mat, y, smooth, bandwidth, n.grid, hc)
    out[i,] <- c(y, bounds$bounds)
  }

  colnames(out) <- c("y", "lower", "upper")
  dte_output <- list("dte" = out)
  class(dte_output) <- c("dte", "list")
  return(dte_output)
}


#' Plot Distribution of Treatment Effects
#'
#' Method for plotting DTE.
#' @param dte_obj Object of class \code{dte} computed by \code{dte_bounds()}.
#' @export
#' @return A \code{ggplot2} plot of the DTE.
#' @examples
#' # Generate data
#' Treat.Mat <- matrix(runif(200) > 0.8, 10, 20)*1
#' Control.Mat <- matrix(runif(100) > 0.8, 10, 10)*1
#'
#' # Generate DTE and plot it
#' bounds <- dte_bounds(Treat.Mat, Control.Mat, -5, 5, smooth = TRUE)
#' plot(bounds)
plot.dte <- function(dte_obj){
  plotdf <- data.frame(dte_obj$dte)
  colnames(plotdf) <- c("y", "lower", "upper")

  plot_obj <- ggplot(plotdf) +
    geom_ribbon(aes(x = y, ymin = lower, ymax = upper), fill="#56B4E9", alpha=0.6, outline.type = "both") + xlab("x") + ylab(expression(paste("P(", DTE <= x,")")))+
    ggtitle("Bounds for DTE") + geom_line(aes(x = y, y = lower), alpha = 0.6)+ geom_line(aes(x = y, y = upper), alpha = 0.6)+
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5), axis.line=element_line())
  print(plot_obj)
  return(plot_obj)
}


# Compute Spectral Treatment Effects
# Output is an STE object
# Makarov Bounds for Symmetric and Bipartite Matrices
#' Spectral Treatment Effects
#'
#' This function computes spectral treatment effects as described in Section 3.2.2 of Auerbach and Cai (2022).
#' @param Treat.Mat Matrix of outcomes for treated units
#' @param Control.Mat Matrix of outcomes for control units
#' @param hc Logical indicating whether or not to incorporate row and column heterogeneity. See Section 5.2 of Auerbach and Cai (2022).
#' @return An object of class \code{ste} containing \tabular{ll}{
#'    \code{stt} \tab STE on treated units sorted in increasing order.\cr
#'    \tab \cr
#'    \code{stc} \tab STE on control units sorted in increasing order.\cr
#'    \tab \cr
#'    \code{treat.cf} \tab Counterfactual for treated units based on eigenvalue imputation. \cr
#'    \tab \cr
#'    \code{control.cf} \tab Counterfactual for control units based on eigenvalue imputation. \cr
#' }
#' @export
#' @examples
#' # Generate Data
#' temp <- matrix(rnorm(100), 10,10)
#' temp <- (temp + t(temp))/2
#' Treat.Mat <- temp*1
#' diag(Treat.Mat) <- 0
#'
#' temp <- matrix(rnorm(400), 20,20)
#' temp <- (temp + t(temp))/2
#' Control.Mat <- temp*1
#' diag(Control.Mat) <- 0
#' # Compute STE
#' effects <- ste(Treat.Mat, Control.Mat)
#' # Compute CDF of STE and plot them
#' effects_cdf <- ste_cdf(effects)
#' plot(effects_cdf$stt)
#' plot(effects_cdf$stc)
ste <- function(Treat.Mat, Control.Mat, hc = FALSE){

  bipartite_flag <- FALSE
  if(!isSquare(Treat.Mat) || !isSquare(Control.Mat)){
    print("Input non-square. They will be treated as bipartite graphs.")
    bipartite_flag <- TRUE
  } else if(!isSymmetric(Treat.Mat) || !isSymmetric(Control.Mat)){
    print("Input asymmetric. They will be treated as bipartite graphs.")
    bipartite_flag <- TRUE
  }

  if(bipartite_flag){
    blowup <- joint_blowup(Treat.Mat, Control.Mat)
    blowup.dim <- dim(blowup[[1]])
    Treat.Mat <- blowup[[1]]
    Control.Mat <- blowup[[2]]
  }

  if(hc){
    Treat.demeaned <- matrix_demean(Treat.Mat)
    Control.demeaned <- matrix_demean(Control.Mat)
    treat.use <- matrix_symmetrize(Treat.demeaned$demeaned.mat, force = bipartite_flag)
    control.use <- matrix_symmetrize(Control.demeaned$demeaned.mat, force = bipartite_flag)
  } else {
    treat.use <- matrix_symmetrize(Treat.Mat, force = bipartite_flag)
    control.use <- matrix_symmetrize(Control.Mat, force = bipartite_flag)
  }

  n.treat <- dim(treat.use)[1]
  n.control <- dim(control.use)[1]
  eig.treat <- eigen(treat.use)
  eig.control <- eigen(control.use)

  # Imputation method doesn't matter for bipartite graphs because of blow-up
  continuous.cf <- joint_blowup_eigen(eig.treat$values, eig.control$values, method = "continuous")

  # Construct Counterfactual by inputing eigenvalues
  # If hc, then this part is from the demeaned portion only
  # Counterfactual of the symmetric embedding
  treat.cf <- eig.treat$vector%*%diag(continuous.cf[[1]])%*%t(eig.treat$vector)
  control.cf <- eig.control$vector%*%diag(continuous.cf[[2]])%*%t(eig.control$vector)

  if(bipartite_flag){
    treat.cf <- treat.cf[1:blowup.dim[1], dim(treat.cf)[1] - blowup.dim[2]:1 + 1]
    control.cf <- control.cf[1:blowup.dim[1], dim(control.cf)[1] - blowup.dim[2]:1 + 1]
  }

  if(hc){

    # Counterfactual fixed effects
    rowfe.cf <- fe.cf(Treat.demeaned$row.means, Control.demeaned$row.means)
    colfe.cf <- fe.cf(Treat.demeaned$col.means, Control.demeaned$col.means)

    rowfe.mat <- list(matrix(Treat.demeaned$row.means, dim(Treat.Mat)[1], dim(Treat.Mat)[2], byrow = F),
                      matrix(Control.demeaned$row.means, dim(Control.Mat)[1], dim(Control.Mat)[2], byrow = F))

    colfe.mat <- list(matrix(Treat.demeaned$col.means, dim(Treat.Mat)[1], dim(Treat.Mat)[2], byrow = T),
                      matrix(Control.demeaned$col.means, dim(Control.Mat)[1], dim(Control.Mat)[2], byrow = T))

    rowfe.cf.mat <- list(matrix(rowfe.cf[[1]], dim(Treat.Mat)[1], dim(Treat.Mat)[2], byrow = F),
                         matrix(rowfe.cf[[2]], dim(Control.Mat)[1], dim(Control.Mat)[2], byrow = F))

    colfe.cf.mat <- list(matrix(colfe.cf[[1]], dim(Treat.Mat)[1], dim(Treat.Mat)[2], byrow = T),
                         matrix(colfe.cf[[2]], dim(Control.Mat)[1], dim(Control.Mat)[2], byrow = T))

    treat.cf.final <- treat.cf + rowfe.cf.mat[[1]] + colfe.cf.mat[[1]]
    control.cf.final <- control.cf + rowfe.cf.mat[[2]] + colfe.cf.mat[[2]]

    stt <- Treat.Mat - treat.cf.final
    stc <- Control.Mat - control.cf.final

  } else {

    treat.cf.final <- treat.cf
    control.cf.final <- control.cf

    stt <- Treat.Mat - treat.cf.final
    stc <- Control.Mat - control.cf.final
  }

  if(bipartite_flag){
    stt <- sort(stt, decreasing = FALSE)
    stc <- sort(stc, decreasing = FALSE)
  } else {
    stt <- sort(stt[upper.tri(stt)], decreasing = FALSE)
    stc <- sort(stc[upper.tri(stc)], decreasing = FALSE)
  }

  out <- list("stt" = stt,
              "stc" = stc,
              "treat.cf" = treat.cf.final,
              "control.cf" = control.cf.final)
  class(out) <- "ste"
  return(out)

}


#' Construct CDFs of Spectral Treatment Effects
#'
#' This function computes CDFs of spectral treatment effects computed by \code{ste()}. The output of this function can be plotted directly using \code{plot()}.
#' @param ste_obj Object of clas 'ste' computed by the function \code{ste()}.
#' @param smooth Logical indicating whether the smoothing should be performed. Smoothing is recommended if potential outcomes are estimated.
#' @param bandwidth 2x1 Array indicating bandwidth to be used for smoothing STT and STC respectively.
#' @param start.ext Number of standard deviations by which to extend the start point of the CDF
#' @param end.ext Number of standard deviations by which to extend the end point of the CDF
#' @param n.grid 2x1 Array indicating the number of grids used to compute the CDF for treatment and control respectively.
#' @return A list \code{ste} containing \tabular{ll}{
#'    \code{stt_cdf} \tab An \code{ste_cdf} object containing the CDF of the STT \cr
#'    \tab \cr
#'    \code{stc_cdf} \tab An \code{ste_cdf} object containing the CDF of the STC \cr
#' }
#' @export
#' @examples
#' # Generate Data
#' temp <- matrix(rnorm(100), 10,10)
#' temp <- (temp + t(temp))/2
#' Treat.Mat <- temp*1
#' diag(Treat.Mat) <- 0
#'
#' temp <- matrix(rnorm(400), 20,20)
#' temp <- (temp + t(temp))/2
#' Control.Mat <- temp*1
#' diag(Control.Mat) <- 0
#' # Compute STE
#' effects <- ste(Treat.Mat, Control.Mat)
#' # Compute CDF of STE and plot them
#' effects_cdf <- ste_cdf(effects)
#' plot(effects_cdf$stt)
#' plot(effects_cdf$stc)
ste_cdf <- function(ste_obj, smooth = FALSE, bandwidth = NULL, start.ext = 1, end.ext = 1, n.grid = c(NULL, NULL)){

  if(class(ste_obj) != "ste"){
    stop("This function can only be used with 'ste' objects.")
  }

  if(smooth & is.null(bandwidth)){
    bandwidth_stt <- (length(ste_obj$stt))^(-1/5)
    bandwidth_stc <- (length(ste_obj$stc))^(-1/5)
  } else if (!smooth){
    bandwidth_stt <- 0
    bandwidth_stc <- 0
  } else {
    bandwidth_stt <- bandwidth[1]
    bandwidth_stc <- bandwidth[2]
  }

  stt_cdf <- smooth_cdf(ste_obj$stt, bandwidth_stt, start.ext, end.ext, n.grid[1])
  class(stt_cdf) <- c("ste_cdf", "array")
  stc_cdf <- smooth_cdf(ste_obj$stc, bandwidth_stc, start.ext, end.ext, n.grid[2])
  class(stc_cdf) <- c("ste_cdf", "array")


  out <- list("stt" = stt_cdf, "stc" = stc_cdf)

  return(out)

}

smooth_cdf <- function(values, bandwidth, start.ext = 1, end.ext = 1, n.grid = NULL){

  if(is.null(n.grid)){
    n.grid <- length(values)
  }

  start <- min(values)-start.ext*sd(values)
  end <- max(values)+end.ext*sd(values)

  grid <-seq(start, end, length.out = n.grid)
  cdf <- grid*0

  if(bandwidth > 0){
    for(i in 1:length(grid)){
      cdf[i] <- mean(ind_smooth(values - grid[i], bandwidth))
    }
  } else {
    for(i in 1:length(grid)){
      cdf[i] <- mean(values - grid[i] <= 0)
    }
  }

  out <- cbind(grid, cdf)
  colnames(out) <- c("y", "cdf.y")
  return(out)
}

#' Plot CDFs of Spectral of Treatment Effects
#'
#' This function plots CDFs of spectral treatment effects computed by \code{ste_cdf()}.
#' @param ste_cdf_obj Object of class \code{ste_cdf} computed by the function \code{ste_cdf()}
#' @return A \code{ggplot2} plot of the STE.
#' @export
#' @examples
#' # Generate Data
#' temp <- matrix(rnorm(100), 10,10)
#' temp <- (temp + t(temp))/2
#' Treat.Mat <- temp*1
#' diag(Treat.Mat) <- 0
#'
#' temp <- matrix(rnorm(400), 20,20)
#' temp <- (temp + t(temp))/2
#' Control.Mat <- temp*1
#' diag(Control.Mat) <- 0
#' # Compute STE
#' effects <- ste(Treat.Mat, Control.Mat)
#' # Compute CDF of STE and plot them
#' effects_cdf <- ste_cdf(effects)
#' plot(effects_cdf$stt)
#' plot(effects_cdf$stc)
plot.ste_cdf <- function(ste_cdf_obj){
  plotdf <- data.frame(ste_cdf_obj)
  plot_obj <- ggplot(plotdf, aes(x = y, y = cdf.y)) + geom_line() +
    xlab("x") + ylab(expression(paste("P(", STE <= x,")"))) +
    ggtitle("Distribution of STE") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5), axis.line=element_line())
  print(plot_obj)
  return(plot_obj)
}

