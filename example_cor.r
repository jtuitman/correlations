######################## EXAMPLES cor.r #########################

n <- 10       # number of variables
rho <- 0.4    # average correlation, 

# Example 1: forward uniform example

l1 <- gen_li_forward_uniform(n,rho)
C1 <- gen_C_from_li(n,l1)
#print(C1)
det(C1)^(1/n)
test_cor(n,rho,C1)

# Example 2: backward uniform example

l2 <- gen_li_backward_uniform(n,rho)
C2 <- gen_C_from_li(n,l2)
#print(C2)
det(C2)^(1/n)
test_cor(n,rho,C2)

# Example 3: forward gaussian example

l3 <- gen_li_forward_gaussian(n,rho)
C3 <- gen_C_from_li(n,l3)
#print (C3)
det(C3)^(1/n)
test_cor(n,rho,C3)

# Example 4: backward gaussian example

l4 <- gen_li_backward_gaussian(n,rho)
C4 <- gen_C_from_li(n,l4)
#print (C4)
det(C4)^(1/n)
test_cor(n,rho,C4)

