#builds the matrices and vectors required for the quadratic form
give_quad <- function(H,j)
{
  A <- H[-j,-j] #the matrix A in quadratic form
  b <- H[-j,j] #the vector b
  list(A, b)
}

#--------- nodewise Lasso function ---------------
# Gives a vector of the inverse hessian matrix

#library(olpsR)

Gr <- function(Sx,Sy,Sxy,Syx,a,b,g1,g2)
{
  l1 <- Sxy%*%b+ 2*g1*Sx%*%a # Del L(alpha)
  l2 <- Syx%*%a +2*g2*Sy%*%b #Del L(beta)
  l3 <- sum(a*Sx%*%a)-1 # Del L(gamma)
  l4 <- sum(b*Sy%*%b)-1
  rbind(l1, l2, l3, l4)
}

nodewise <- function(j, H, lam)
{
  l <- give_quad(H, j)
  A <- l[[1]]
  b <- l[[2]]

  init_val <- opt_init(l, lam)
  init_val <- pp(init_val)

  #making A-tilde
  t1 <- cbind(A, -A)
  t2 <- cbind(-A, A)
  tA <- rbind(t1, t2)
  #making b-tilde
  tb <- rbind(2*b-lam, -2*b-lam)

  B <-  10/lam
  # Writing the bounds as linear constraint
  l <- length(init_val)
  conm <- diag(rep(1, l))
  conm <- rbind(conm, -conm)
  b <- c(rep(-10^(-10), l), rep(-B, l))

  res <- stats::constrOptim(init_val, f=nodewise_obj, grad=nodewise_grad, A=tA, b=tb,
                method="CG", ui=conm, ci=b)

  sol <- res$par
  conv <- res$convergence
  # Projecting on a simplex to ensure that the L-1 norm of the solution
  #is bounded above
  #by b
  #sol <- projsplx(sol, b = 10/lam)

  #Going back to the original problem
  m <- length(sol)
  sol1 <- as.vector(sol[1 : (m/2)])
  sol2 <- as.vector(sol[(m/2+1) : m])
  sol <- sol1 - sol2
  sol <- matrix(sol,ncol=1)

  if(j==1) {Tau <- c(1,-sol)} else {
    Tau <- c(-sol[1:(j-1)], 1, -sol[j : length(sol)])
  if(j==(length(sol)+1)) Tau <- c(-sol, 1)

  }

  tau <- sum(as.vector(H%*%Tau)*Tau)+lam*sum(abs(sol))/2

  #setting small numbers to 0
  v <- Tau/tau
  v[v<10^(-10)] <- 0
  return(v)

}

#-------- objective function of optim --------------

nodewise_obj <- function(x, A, b)
{
  t1 <- A%*%x
  sum(as.vector(t1*x)- as.vector(b*x))
}

#---------- gradient function for optim -------------

nodewise_grad <- function(x, A, b)
{

 as.vector( 2*A%*%x)-as.vector(b)
}

#------- gives positive part and negative part ------

pp <- function(x)
{
  t1 <- ifelse(x>0, x, 0)
  t2 <- ifelse(x<0, -x, 0)
  c(t1,t2)
}
