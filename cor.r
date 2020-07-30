gen_li_forward_uniform <- function(n,rho){
  
  L <- sqrt(n+rho*(n^2-n))
  
  li <- 1
  if(2 <= n-2){
    for (i in 2:(n-2)){
      a <- max(L-(n-i),abs(li[i-1]-1))
      b <- min(L+(n-i),li[i-1]+1)
      li[i] <- runif(1,min=a,max=b) # uniform distribution
    } 
  } 
  a <- max(abs(li[n-2]-1),abs(L-1))
  b <- min(li[n-2]+1,L+1)
  li[n-1] <- runif(1,min=a,max=b)
  li[n] <- L
  return(li)
}


gen_li_backward_uniform <- function(n,rho){

  li <- 0
  for (i in 2:n){
    li[i] <- 0 # initialisation
  }

  L <- sqrt(n+rho*(n^2-n))

  li[n] <- L
  if(2 <= n-1){
    for (i in (n-1):2){
      a <- abs(li[i+1]-1)
      b <- min(li[i+1]+1,i)
      li[i] <- runif(1,min=a,max=b) 
    }
  }
  li[1] <- 1
  return(li)
}


gen_li_forward_gaussian <- function(n,rho){

  L <- sqrt(n+rho*(n^2-n))

  li <- 1
  if(2 <= n-2){
    for (i in 2:(n-2)){
      a <- max(L-(n-i),abs(li[i-1]-1))
      b <- min(L+(n-i),li[i-1]+1)
      ex <- li[i-1]+(L-li[i-1])/(n-i+1)
      if(ex >= a & ex <= b){
        done <- FALSE
        while(done == FALSE){
          r=rnorm(1,ex,min(ex-a,b-ex)/2)
          if(r >= a & r <= b){
            done <- TRUE
          }
        }
    } else{
        done <- FALSE
        while(done == FALSE){
          r=rnorm(1,(a+b)/2,(b-a)/4)
          if(r >= a & r <= b){
            done <- TRUE
          }
        }
      }
      li[i] <- r
    }     
  }     
  a <- max(abs(li[n-2]-1),abs(L-1))
  b <- min(li[n-2]+1,L+1)
  done <- FALSE
  while(done == FALSE){
    r=rnorm(1,(a+b)/2,(b-a)/4)
    if(r >= a & r <= b){
      done <- TRUE
    }
  } 
  li[n-1] <- r
  li[n] <- L
  return(li)
}


gen_li_backward_gaussian <- function(n,rho){

  li <- 0
  for(i in (1:n)){
    li[i] <- 0 # initialisation
  }
  
  L <- sqrt(n+rho*(n^2-n))
  
  li[n] <- L
  
  if(2 <= n-1){
    for(i in ((n-1):2)){
      a <- abs(li[i+1]-1)
      b <- min(li[i+1]+1,i)
      ex <- li[i+1]-(li[i+1]-1)/i
      if(ex >=a & ex <= b){
        done <- FALSE
        while(done == FALSE){
          r=rnorm(1,ex,min(ex-a,b-ex)/2)
          if(r >= a & r <= b){
            done <- TRUE
          }
        }
      }
      else{
        done <- FALSE
        while(done == FALSE){
          r=rnorm(1,(a+b)/2,(b-a)/4)
          if(r >= a & r <= b){
            done <- TRUE
          }
        }   
      }
      li[i] <- r
    }
  } 
  
  li[1] <- 1

  return(li)

}


gen_T_from_li <- function(n,li){
  
  rows <- vector(mode="list", length=n)
  s <- vector(mode="list", length=n)

  x <- rnorm(n)
  x <- x/sqrt(sum(x^2)) # random vector of length 1
  rows[[1]] <- x 
  s[[1]] <- x

  for (i in 2:n){
    p <- (li[i]^2-li[i-1]^2-1)/2
    x <- rnorm(n)
    y <- x - s[[i-1]]*sum(x*s[[i-1]])/sum(s[[i-1]]^2) # y is a random vector orthogonal to s[[i-1]]
    z <- s[[i-1]]*p/sum(s[[i-1]]^2)                   # <z,s[[i-1]]> is already what it should be, 
                                                      # but z not of length 1
    if(sum(z^2)<=1){
      ti <- z+y*sqrt((1-sum(z^2))/sum(y^2))           # add the right multiple of y to z, so that 
                                                      # ti is of length 1
    }
    else{
      ti <- z
    }
    rows[[i]]<-ti
    s[[i]]<-s[[i-1]]+ti
  }  
  
  T <- matrix(,n,n)
  for (i in 1:n){
    for (j in 1:n){
      T[i,j] <- rows[[i]][j] # the matrix T with rows given by 'rows' 
    }
  }  

  return(T)
}


gen_C_from_T <- function(T){
  
  C <- T %*% t(T) # correlation matrix C=T*transpose(T)
   
  return(C) 

}


gen_C_from_li <- function(n,li){

  T <- gen_T_from_li(n,li)
  C <- T %*% t(T) # correlation matrix C=T*transpose(T)
  
  return(C)

}


gen_C_from_rho <- function(n,rho){

  li <- gen_li_forward_gaussian(n,rho)
  C <- gen_C_from_li(n,li)
  
  return(C)

}


test_cor <- function(n,rho,C){

  S <- 0
  for (i in 1:n){
    for (j in 1:n){
      if(i != j){S <- S + C[i,j]}
    }
  }
  return(S-rho*(n^2-n)) # should be close to zero
}

