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