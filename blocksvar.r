library(far)


my_sum <- function(v,i,j){

  s <- 0
  
  if(i<=j){
    for (k in i:j){
      s <- s + v[k]
    }
  }

  return(s)
}


gen_li_forward_gaussian <- function(n,sigmai,S){


  li <- sigmai[1]
  if(2 <= n-1){
    for (i in 2:(n-1)){
      a <- max(S-my_sum(sigmai,i+1,n),sigmai[n]-my_sum(sigmai,i+1,n-1)-S,abs(li[i-1]-sigmai[i]))
      b <- min(S+my_sum(sigmai,i+1,n),li[i-1]+sigmai[i])
      ex <- li[i-1]+sigmai[i]*(S-li[i-1])/(my_sum(sigmai,i,n))
      if(ex >= a & ex <= b){
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
    }
    else{
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


random_On <- function(n){
 
  A=matrix(,n,n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j]=rnorm(1)
    }
  }
  
  qrdata=qr(A)
  Q=qr.Q(qrdata)
  R=qr.R(qrdata)
 
  for(i in 1:n){
    if(R[i,i] < 1){
      for(j in 1:n){
        Q[j,i]=-Q[j,i]
      }
    }
  }
  
  return(Q)
 
}


random_On_v_to_w <- function(n,v,w){
  
  V = orthonormalization(v) # Gramm-Schmidt orthogonalization from package "far"
  W = orthonormalization(w)
  
  g = random_On(n-1) 
  
  M = matrix(0,nrow=n,ncol=n) # initialisation
  M[1,1] = 1
  for(i in 2:n){
    for(j in 2:n){
      M[i,j] = g[i-1,j-1]
    }
  }
  
  return(W %*% M %*% solve(V)) # solve is matrix inverse, confusing R syntax!

}


gen_w_from_Si_S <-function(k,Si,S){

  li = gen_li_forward_gaussian(k,Si,S)
  
  
  W = gen_T_from_li(k,Si,li)
  
  for(i in 1:k){
    for(j in 1:k){
      W[i,j]=Si[i]*W[i,j]
    }
  }
  
  return(W)
}


gen_C_from_Ti_S<-function(Ti,sigmaij,S){

  k=length(Ti)
  n <- 0
  for(i in 1:k){
    n = n + nrow(Ti[[i]])
  }

  V=matrix(0,nrow=k,ncol=n)
  Si=0
  for(i in 1:k){
    vec=0
    for(j in 1:nrow(Ti[[i]])){
      vec=vec+sigmaij[[i]][j]*Ti[[i]][j,]
    }
    for(j in 1:nrow(Ti[[i]])){
      V[i,j]=vec[j]
    }
    Si[[i]]=sqrt(sum(vec^2))
  }
  
  W = gen_w_from_Si_S(k,Si,S)
  
  Wnew=matrix(0,nrow=k,ncol=n)
  for(i in 1:k){
    for(j in 1:k){
      Wnew[i,j] = W[i,j]
    }
  }
  
  M = random_On(n)
  W = Wnew %*% M
  
  T = matrix(0,nrow=n,ncol=n)
  ct=0
  for(i in 1:k){
    Mi=random_On_v_to_w(n,V[i,],W[i,])
    for(j in 1:nrow(Ti[[i]])){ 
      ct = ct+1
      tij=matrix(0,nrow=n,ncol=1)
      for(l in 1:nrow(Ti[[i]])){
        tij[l,1]=Ti[[i]][j,l]
      }
      T[ct,]=Mi%*%tij     
    }
  }
  
  C = T%*%t(T)
  
  return(C)

}


gen_C_from_Ci_S<-function(Ci,sigmaij,S){

  k = length(Ci)
  Ti = Ci

  for(i in 1:k){
    Ti[[i]] = t(chol(Ci[[i]]))
  }  
  
  C = gen_C_from_Ti_S(Ti,sigmaij,S)

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
