my_sum <- function(v,i,j){

  s <- 0
  
  if(i<=j){
    for (k in i:j){
      s <- s + v[k]
    }
  }

  return(s)
}


gen_li_forward_uniform <- function(n,sigmai,S){
  
  li <- sigmai[1]
  if(2 <= n-1){
    for (i in 2:(n-1)){
      a <- max(S-my_sum(sigmai,i+1,n),sigmai[n]-my_sum(sigmai,i+1,n-1)-S,abs(li[i-1]-sigmai[i]))
      b <- min(S+my_sum(sigmai,i+1,n),li[i-1]+sigmai[i])
      li[i] <- runif(1,min=a,max=b) 
    }
  } 
  li[n] <- S
  return(li)
}


gen_li_backward_uniform <- function(n,sigmai,S){

  li <- 0
  for (i in 2:n){
    li[i] <- 0 # initialisation
  }

  li[n] <- S
  if(2 <= n-1){
    for (i in (n-1):2){
      a <- max(sigmai[1]-my_sum(sigmai,2,i),abs(li[i+1]-sigmai[i+1]))
      b <- min(li[i+1]+sigmai[i+1],my_sum(sigmai,1,i))
      li[i] <- runif(1,min=a,max=b) 
    }
  } 
  li[1] <- sigmai[1]
  return(li)
}


gen_li_forward_gaussian <- function(n,sigmai,S){

  li <- sigmai[1]
  
  if(2 <= n-1){
    for (i in 2:(n-1)){
      a <- max(S-my_sum(sigmai,i+1,n),sigmai[n]-my_sum(sigmai,i+1,n-1)-S,abs(li[i-1]-sigmai[i]))
      b <- min(S+my_sum(sigmai,i+1,n),li[i-1]+sigmai[i])
      ex <- li[i-1]+sigmai[i]*(S-li[i-1])/(my_sum(sigmai,i,n))
      if(ex >= a & ex <= b){i
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
  li[n] <- S
  return(li)
}


gen_li_backward_gaussian <- function(n,sigmai,S){

  li <- 0
  for(i in (1:n)){
    li[i] <- 0 # initialisation
  }
  
  li[n] <- S
  
  if(2 <= n-1){
    for(i in ((n-1):2)){
      a <- max(abs(li[i+1]-sigmai[i+1]),sigmai[1]-my_sum(sigmai,2,i))  
      b <- min(li[i+1]+sigmai[i+1],my_sum(sigmai,1,i))
      ex <- li[i+1]-(li[i+1]-sigmai[1])/my_sum(sigmai,1,i)
      if(ex >=a & ex <= b){
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
  li[1] <- sigmai[1]

  return(li)

}


gen_T_from_li <- function(n,sigmai,li){

  rows <- vector(mode="list", length=n)
  s <- vector(mode="list", length=n)

  x <- rnorm(n)
  x <- x/sqrt(sum(x^2)) # random vector of length 1
  rows[[1]] <- x 
  s[[1]] <- sigmai[1]*x
 
  
  for (i in 2:n){
    p <- (li[i]^2-li[i-1]^2-sigmai[i]^2)/(2*sigmai[i])
    x <- rnorm(n)
    y <- x - s[[i-1]]*sum(x*s[[i-1]])/sum(s[[i-1]]^2) # y is a random vector orthogonal to s[[i-1]]
    z <- s[[i-1]]*p/sum(s[[i-1]]^2)                   # <z,s[[i-1]]> is already what it should be, 
                                                      # but z not of length 1
    if(sum(z^2)<=1){
      ti <- z+y*sqrt((1-sum(z^2))/sum(y^2))           # add the right multiple of y to z, so that 
                                                      # ti is of length 1
    } else{
      ti <- z
    }
    rows[[i]]<-ti
    s[[i]]<-s[[i-1]]+sigmai[i]*ti
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


gen_C_from_li <- function(n,sigmai,li){

  T <- gen_T_from_li(n,sigmai,li)
  C <- T %*% t(T) # correlation matrix C=T*transpose(T)
  
  return(C)
  
}


gen_C_from_S <- function(n,sigmai,S){

  li <- gen_li_forward_gaussian(n,sigmai,S)
  C <- gen_C_from_li(n,sigmai,li)
  
  return(C)

}


test_cor <- function(n,sigmai,C,S){

  s <- 0
  for (i in 1:n){
    for (j in 1:n){
      s <- s+C[i,j]*sigmai[i]*sigmai[j]
    }
  }
  return(S-abs(s)^(1/2)) # should be close to zero
}



