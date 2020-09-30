# find the numerical root of g(x)=0 on [a,b]
# also return the number of bisections
root_bis <- function(g,a,b,eps=10^(-6))
{
  # check if a root exists on [a,b]
  if (abs(g(a))<eps) {return(list('root'=a,'num_bisect'=0))} else
    if (abs(g(b))<eps) {return(list('root'=b,'num_bisect'=0))} else
      if (g(a)*g(b)>0) {return(paste('No root on [',a,',',b,']'))}
  
  # find root on (a,b) by bisection
  left <- a
  right <- b
  count <- 0
  while (right-left >= eps)
  {
    mid <- (left + right) / 2
    if (g(left)*g(mid)<0) {right <- mid} else {left <- mid}
    count <- count + 1
  }
  # output numerical root on [a,b] and the number of bisections
  return(list("root"=mid, "num_bisect"=count))
}



# find numerical root of g(x)=0 starting from initial guess
# also return the number of Newton-Raphson iterations
root_NR <- function(g, init, eps=10^(-6))
{
	xx <- init
	count <- 0
	while(abs(g(xx)) >= eps)
	{
		xx <- xx - g(xx) * 2*eps/(g(xx+eps) - g(xx-eps))
		count <- count + 1
	}
	return(list('root'=xx, 'num_iter'=count))
}

# fit X to a Gaussian Mixture Model defined by
# component weights, component means, and component sd

EM_init <- function(X,w_init,mu_init,sigma_init,eps=10^(-6))   # initialize by prior knowledge
{
  theta_old <- 0
  theta_new <- c(w_init, mu_init, sigma_init)
  
  while (sum((theta_new-theta_old)^2) >= eps)
  {
    theta_old <- theta_new
    # construct matrix Omega
    T <- matrix(theta_old, ncol=3)
    OMG <- c()
    for (i in 1:length(X))
    {
      omg <- T[,1]*dnorm(X[i], mean=T[,2], sd=T[,3])
      OMG <- c(OMG, omg/sum(omg))
    }
    OMG <- t(matrix(OMG, ncol=length(X)))
    # update
    w <- colSums(OMG)
    w <- w/sum(w)
    mu <- as.vector(X%*%OMG / colSums(OMG))
    sig <- diag(t(outer(X,mu,`-`)^2) %*% OMG) / colSums(OMG)
    theta_new <- c(w,mu,sqrt(as.vector(sig)))
  }
  
  # output parameters as a matrix
  return(matrix(theta_new, ncol=3))
}

EM_diff <- function(X,K)  # initialize by overall averge and sd of sample
{
  m <- mean(X)
  s <- sd(X)
  EM_init(X,
          w_init = rep(1/K, K),
          mu_init = seq(m-0.5*s*(K-1), m+0.5*s*(K-1), length.out=K),
          sigma_init = rep(s,K)
  )
}

EM_group <- function(X,K) # intialize by K-means group
{
  clusters <- kmeans(X,K)$cluster
  EM_init(X, 
          w_init = table(clusters)/length(X),
          mu_init = tapply(X, clusters, mean),
          sigma_init = tapply(X, clusters, sd)
  )
}