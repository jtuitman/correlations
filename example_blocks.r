######################## EXAMPLES blocks.r #########################

C1 = matrix(c(1.0000000, 0.3581445, 0.6069251, 0.3581445, 1.0000000, 0.5349304, 0.6069251, 0.5349304, 1.0000000),nrow=3,ncol=3)
C2 = matrix(c(1.0000000, 0.8115993, 0.1405289, 0.8115993, 1.0000000, 0.0478718, 0.1405289, 0.0478718, 1.0000000),nrow=3,ncol=3)
C3 = matrix(c(1.00000000, -0.01722136, 0.05192323, -0.01722136,  1.00000000, 0.71529813, 0.05192323,  0.71529813, 1.00000000),nrow=3,ncol=3)

# C1,C2 and C3 are 3x3 correlation matrices with average correlation 1/2,1/3,1/4, respectively.

print(C1)
print(C2)
print(C3)

Ci = list(C1,C2,C3)

C = gen_C_from_Ci(Ci) # 9x9 correlation matrix with blocks C1,C2,C3 along the diagonal

print(C)

