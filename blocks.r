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


gen_C_from_Ti <- function(Ti){

  k=length(Ti)
  n <- 0
  for(i in 1:k){
    n = n + nrow(Ti[[i]])
  }
  
  T<-matrix(,n,n)
  cnt<-0
  
  for(i in 1:k){
    Mi <- random_On(n)
    Tii <- Ti[[i]]
    ni <- nrow(Tii)
    for(j in 1:ni){
      cnt=cnt+1
      tij=Tii[j,]
      for(l in (ni+1):n){
        tij[l]=0
      }
      for(l in 1:n){
        T[cnt,l]=(Mi%*%tij)[l]
      } 
      
    }
  }

  C <- T%*%t(T)
  
  return(C)
  
}


gen_C_from_Ci <- function(Ci){

  k=length(Ci)
  Ti=Ci

  for(i in 1:k){
    Ti[[i]] = t(chol(Ci[[i]]))
  }  
  
  C = gen_C_from_Ti(Ti)
  
  return(C)
}

