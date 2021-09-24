# Calculating Hessian
# Input: Sxy = Cross-covariance matrix
## Sx, Sy = Marginal covariance matrices
## x, y = estimates of x and y vectors
## rho_init = initial estimate of the first canoncial correlation
## C = Hessian parameter

 hessian <- function(Sxy, Sx, Sy, x, y, rho_init, C)
 {
   Syx <- t(Sxy)
   H1 <- 2*rho_init*Sx + 2*C*(Sx %*% (x %o% x) %*% Sx)

   H12 <- 2*(2-C) * (Sx %*% (x %o% y) %*% Sy) - 2*Sxy

   H2 <- 2*rho_init*Sy + 2*C*(Sy %*% (y %o% y) %*% Sy)

   H21 <- t(H12)

   temp1 <- cbind(H1, H12)
   temp2 <- cbind(H21, H2)
   rbind(temp1, temp2)
 }




