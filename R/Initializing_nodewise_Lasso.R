#Nodewise Lasso initialization

# input matrix A is the Hessian.
#--------------- Approximating A by a positive definite matrix ---------------------
give_approx <- function(A)
{
  egn <- eigen(A, symmetric=TRUE)
  V <- egn$vectors
  Lam <- egn$values

  if(min(Lam)<0)
    Lam <- Lam-min(Lam) + 0.3

  apA <- V%*%diag(Lam)%*%solve(V)
 Matrix::forceSymmetric(apA)
}

#-----------------------------------------------------------------------------

# gives an initial solution
# The input l is a list of a matrix A (non p.d. but symmetric) and a vector b
# lam is the lambda, sqrt(log p/n)
opt_init <- function(l,lam)
{
  A <- give_approx(l[[1]])
  b <- l[[2]]
  x <- Variable(length(b))

  obj <- sum(quad_form(x,A),-t(b)%*%x,lam*abs(x))
  prob <- Problem(Minimize(obj))

  # Here we do not have the constraint because we have a convex program
  result <- solve(prob)

  val <- result$value # Optimal objective
  gamma <- result$getValue(x) # Optimal variables
  gamma
}




