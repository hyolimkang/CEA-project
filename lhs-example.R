# Latin hypercube sample - example
# https://cran.r-project.org/web/packages/lhs/vignettes/lhs_basics.html

# library (lhs)

# set the seed for reproducibility
set.seed (1)

runs <- 100

# a design with n samples from k parameters
A <- randomLHS (n = runs, 
                k = 4) 
A

# It is common to transform the margins of the design (the columns) 
# into other distributions
lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A)) # 
lhs_sample [,1] <- qnorm (p = A[,1], mean = 0, sd = 1)
lhs_sample [,2] <- qlnorm (p = A[,2], meanlog = 0.5, sdlog = 1)
lhs_sample[,3] <- A[,3]
# B[,3] <- qunif (p = A[,3], min = 0, max = 1)
lhs_sample [,4] <- qunif (p = A[,4], min = 7, max = 10)