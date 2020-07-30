# Example 1: forward uniform example

sigmai <- 1
for(i in 2:9){
  sigmai[i] = 1
} 
sigmai[10] = 3
n=10
S=2


li <- gen_li_forward_uniform(n,sigmai,S)
C <- gen_C_from_li(n,sigmai,li)
test_cor(n,sigmai,C,S)

# Example 2: backward uniform example

sigmai <- 3 
for(i in 2:10){
sigmai[i] = 1
}
n=10
S=0

li <- gen_li_backward_uniform(n,sigmai,S)
C <- gen_C_from_li(n,sigmai,li)
test_cor(n,sigmai,C,S)

# Example 3: forward gaussian example

sigmai <- 1
for(i in 2:9){
  sigmai[i] = 1
} 
sigmai[10] = 5 
n=10
S=3 

li <- gen_li_forward_gaussian(n,sigmai,S)
C <- gen_C_from_li(n,sigmai,li)
test_cor(n,sigmai,C,S)

# Example 4: backward gaussian example

sigmai <- 3 
for(i in 2:10){
sigmai[i] = 1
}
n=10
S=0

li <- gen_li_backward_gaussian(n,sigmai,S)
C <- gen_C_from_li(n,sigmai,li)
test_cor(n,sigmai,C,S)


