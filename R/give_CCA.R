#----------------------
# Inputs:
#----------------------
## l = a list with the first item the initial alpha estimate and the second item the initial beta estimate
## C = the C parameter for the Hessian matrix
## X = the first set of variables
## Y = the second set of variables
## elements = which of the (p + q) parameters should the variance be calculated for
## nlC = the nodewise Lasso tuning parameter
#----------------------
# Outputs a list containing:
#----------------------
## matx = a data frame with debiased and initial x vectors in first and second columns, respectively
## maty = analagous data frame as matx but for the y vectors
## cond = the condition number of the Hessian matrix
## var_calc = a (p + q)-length vector with the variance of the parameters specified by the "elements" argument. Non-specified positions return "NA"
## rhosq_db = the debiased estimate of the squared first canonical correlation
## rhosq_var = the variance estimate of the debiased squared first canonical correlation
#' De-biased CCA estimators
#'
#' Suppose X and Y are random vectors with variance \eqn{\Sigma_x} and
#' \eqn{\Sigma_y},  and cross-covariance matrix \eqn{\Sigma_{xy}}.
#'  The leading canonical correlation coefficient  is defined as
#'        \deqn{\rho=\max_{\alpha,\beta}\alpha'\Sigma_{xy}\beta}
#'        where the maximization is over
#'          \eqn{\alpha: \alpha'\Sigma_x\alpha=1} and \eqn{\beta: \beta'\Sigma_y\beta=1}.
#'         The solutions to the above optimization problems are known as
#'         the canonical directions. If the estimators of \eqn{\alpha} and \eqn{\beta}
#'         are sparse, then a first order bias is incurred. Laha et al. (2021)
#'         introduces a bias correction method, which provides  elementwise \eqn{\sqrt{n}}-consistent
#'         estimators of \eqn{\sqrt{\rho}\alpha} and \eqn{\sqrt{\rho}\beta},
#'         provided the preliminary estimators of \eqn{\alpha} and \eqn{\beta} satisfy some \eqn{l_1} and \eqn{l_2}
#'         consistency properties. We also provide a \eqn{\sqrt{n}}-consistent  estimator of
#'         \eqn{\rho^2}. See Laha et al. (2021) for more details.
#'The following function provides de-biased estimators of \eqn{\sqrt{\rho}\alpha},
#'           \eqn{\sqrt{\rho}\beta}, \eqn{\rho^2}, and estimates of their  variances.
#'
#' @param X A matrix with n rows and p columns, the first dataset.
#' @param Y A matrix with n rows and q columns, the second dataset.
#' @param alpha A vector of length p, estimator of \eqn{\alpha}, the  leading canonical vector corresponding to X.
#' @param beta A vector of length q, estimator of \eqn{\beta}, the  leading canonical vector corresponding to Y.
#' @param C An optional constant. Default is 2.
#' @param elements Optional. A vector of integers taking values in the set {1,2,...,p+q}. Corresponds to variance estimates of
#'                 the debiased estimators of  \eqn{\alpha_i} and \eqn{\beta_j}. See details. The default value
#'                   is the vector {1,2,...,p+q}.
#' @param nlC An optional  constant. Must be positive. The default is \eqn{\sqrt{\log(p+q)/n}}.
#'
#' @details \code{elements} The variance estimates corresponding to
#'                \eqn{\alpha_i} (or \eqn{\beta_j}) are calculated
#'                only if i (or p+i) are in \code{elements}. Otherwise, NA will be returned.
#'@details \code{nlC} The node-wise Lasso parameter \eqn{\lambda^{nl}_j} in Laha et al. (2021).
#'              Should be a constant multiple of  \eqn{\sqrt{\log(p+q)/n}}.
#' @details \code{C} The constant C in Lemma 1 of Laha et al. (2021). The value of 2 has been used in this paper.
#' @details \code{nlC}  Penalty parameter for a
#'                \eqn{l_1} penalty term, required in the nodewise Lasso step,
#'                          during inversion of the hessian matrix.
#'
#'@return  A list. Contains the followings:
#'\describe{
#'\item{matx}{ A data frame with de-biased and initial estimators of \eqn{\alpha}
#'                in first and second columns, respectively.}
#'\item{maty}{ A data frame with de-biased and initial estimators of \eqn{\beta}
#'                in first and second columns, respectively.}
#'\item{var_calc}{ A (p + q)-length vector with the variance of the parameters
#'               specified by the \code{elements} argument. Non-specified
#'               positions return "NA".}
#'\item{rhosq_db}{ The debiased estimate of the squared first canonical correlation.}
#'\item{rhosq_var}{ The variance estimate of the debiased squared first canonical correlation.}
#' }
#'
#' @keywords cca
#'
#' @references Laha, N., Huey, N., Coull, B., and Mukherjee, R. (2021). \emph{On Statistical
#' Inference with High Dimensional Sparse CCA}, submitted.
#'
#' @examples
#' library(mvtnorm)
#' library(dplyr)
#' library(CVXR)
#'
#' #Simulate standard normal data matrix: first generate alpha and beta
#' p <- 50; q <- 50; al <- c(rep(1, 10), rep(0, 40));
#' be <- c(rep(0,25), rnorm(25,1))
#'
#' #Normalize alpha and beta
#' al <- al/sqrt(sum(al^2)); be <- be/sqrt(sum(be^2))
#'
#'  #Set n and rho
#' n <- 300; rho <- 0.5
#'
#' #Creating the covariance matrix
#'Sigma_mat <- function(p,q,al,be, rho)
#' {
#'  Sx <- diag(rep(1,p), p, p)
#'  Sy <- diag(rep(1,q), q, q)
#'  Sxy <- tcrossprod(crossprod(rho*Sx, outer(al, be)), Sy)
#'  Syx <- t(Sxy)
#'  rbind(cbind(Sx, Sxy), cbind(Syx, Sy))
#'}
#'truesigma <- Sigma_mat(p,q,al,be, rho)
#'
#'#Finally simulating the data
#'Z <- mvtnorm::rmvnorm(n, sigma = truesigma)
#'x <- Z[,1:p]
#' y <- Z[,(p+1):(p+q)]
#' elements <- 1:p
#' nlC <- log(p+q)/n
#'
#'  # Preliminary estimators: Mai(2019)'s SCCA estimators
#' temp <- cca.mai(x,y)
#' ha <- temp[[1]]
#' hb <- temp[[2]]
#'
#' #Call give_CCA
#' give_CCA(ha, hb, x, y)
#' @export
give_CCA <- function(alpha, beta,  X, Y, C, elements, nlC)
{
  a <- alpha %>% methods::as("sparseMatrix")
  b <- beta %>% methods::as("sparseMatrix")

  Sx <- stats::var(X)
  Sy <- stats::var(Y)
  Sxy <- stats::cov(X,Y)
  Syx <- t(Sxy)
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(X)

  #Setting default variables
  if(missing(C)) C <- 2
  if(missing(elements)) elements <- 1:(p+q)
  if(missing(nlC)) nlC <- sqrt(log(p+q)/n)

  # Getting an initial rho estimate
  rho_init <- abs(sum(a*(Sxy%*%b)))

  xhat <- sqrt(rho_init)*a
  yhat <- sqrt(rho_init)*b
  xmat <- X
  ymat <- Y

  #Hessian
  H <- hessian(Sxy, Sx, Sy, a, b, rho_init, C = C)
  cond <- 1/suppressWarnings(rcond(H))

  #Setting lambda
  lam = nlC * sqrt(log(p+q)/n)

  #Nodewise Lasso
  list.vec <- lapply(1:ncol(H), nodewise, lam=lam, H=H)

  #The estimate of the precision matrix (inverse of Hessian)
  Phi <- matrix(unlist(list.vec), byrow=TRUE, ncol= p+q)

  # Calculate the pseudovariables Z_j: j = 1...n (for variance calculation)
  Z <- rep(NA, n)
  var_calc <- rep(NA, p + q)

  for(i in elements){
    for (j in 1:n) {
      piece1 <- ((rho_init*t(Phi[(1:p), i]) + ((t(Phi[(1:p), i])%*% Sx %*% xhat)%*%t(as.matrix(xhat)))) %*% (outer(xmat[j,],xmat[j,],"*")%*%xhat))
      piece2 <- t(Phi[(1:p), i])%*% xmat[j,] %*% t(ymat[j,])%*%yhat
      piece3 <- ((rho_init*t(Phi[(p+1):(p+q), i]) + ((t(Phi[(p+1):(p+q), i])%*% Sy %*% yhat)%*%t(as.matrix(yhat)))) %*% (outer(ymat[j,], ymat[j,],"*") %*% yhat))
      piece4 <- t(as.matrix(xhat))%*% xmat[j,] %*% t(ymat[j,]) %*% (Phi[(p+1):(p+q), i])
      Z[j] <- piece1 - piece2 + piece3 - piece4
    }
    var_calc[i] <- stats::var(Z)
  }

  #Gradient
  grad1 <- (2*rho_init*(Sx %*% xhat)) - 2*(Sxy %*% yhat)
  grad2 <- (2*rho_init*(Sy %*% yhat)) - 2*(Syx %*% xhat)
  temp <- rbind(grad1,grad2)

  # De-biasing
  de.bias <- rbind(xhat,yhat) - Phi%*%temp
  xhatdebias <- de.bias[1:p]
  yhatdebias <- de.bias[(p+1):(p+q)]

  matx <- data.frame(d.a=xhatdebias, i.a=as.vector(xhat))
  maty <- data.frame(d.b=yhatdebias, i.b=as.vector(yhat))

  xhat <- as.vector(xhat)
  yhat <- as.vector(yhat)

  #Calculating the de-biased squared first canonical correlation
  rhosq_db <- (t(as.matrix(xhat)) %*% Sxy %*% yhatdebias) + (t(xhatdebias) %*% Sxy %*% yhat) - (t(as.matrix(xhat)) %*% Sxy %*% yhat)

  #Calculating the variance of the above estimate
  rhosq_var <- vector(length = n)
  for (i in 1:n) {
      X_i <- xmat[i,]
      Y_i <- ymat[i,]
      rhosq_var[i] <- ((crossprod(xhat, X_i)^2) * (crossprod(yhat, Y_i)^2)) / n
  }

  rhosq_var <- sum(rhosq_var) - ((rho_init)^4)

  #Putting the output list together
  listestim <- list(matx, maty,  var_calc, rhosq_db, rhosq_var)
  return(listestim)
}

#' Leading canonical correlates  based on Mai and Zhang (2019).
#'
#' The code for the general set up
#'  is provided in Qing Mai's website. We include the special rank one
#'case here for completeness because we use this SCCA estimator in an example.
#'
#'@param xmat The X matrix, a matrix with n rows and p columns.
#'@param ymat The Y matrix, a matrix with n rows and q columns.
#'@details The \eqn{l_1} penalty  parameters corresponding to
#'            \eqn{\alpha} and \eqn{\beta} are
#'           taken to be \eqn{\sqrt{(\log p)/n}} and
#'           \eqn{\sqrt{(\log q)/n}}, respectively.
#'
#'@return  A list. Contains the followings:
#'\describe{
#'\item{sa}{Estimator of \eqn{\alpha}, leading left canonical direction.}
#'\item{sb}{Estimator of \eqn{\beta}, leading right canonical direction.}
#'\item{srho}{Estimator of leading canonical correlation.}
#'  }
#'
#'@references Mai, Q. and Zhang, X. (2019). \emph{An iterative penalized least squares approach to sparse
#' canonical correlation analysis}, Biometrics, 75, 734-744.
#' @export
cca.mai <- function(xmat, ymat)
{
  n <- nrow(xmat)
  p <- ncol(xmat)
  q <- ncol(ymat)
  s.cca <- SCCA(xmat, ymat, lambda.alpha=sqrt(log(p)/n), lambda.beta = sqrt(log(q)/n))
  sa <- s.cca$alpha
  sb <-s.cca$beta
  srho <- sum(t(sa)%*%cov(xmat,ymat)%*%sb)
  list(sa, sb, srho)
}
