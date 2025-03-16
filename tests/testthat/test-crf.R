set.seed(0)
n=2; I=1000; N=n*I
x = matrix(rnorm(2*N),nrow=N,ncol=2)
epsilon = rnorm(N)
sqrtCovY = matrix(c(2,1,1,2)/sqrt(5), 2, 2)
for (i in 1:I) epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
y = tanh(x[,1]) + tanh(x[,2]) + epsilon
data <- data.frame(y=y, x=x, id=rep(1:(N/n),each=n))

test_that("crf works for: Training MSE; equicorr; L=NULL", {
  tst <- crf(y~x.1+x.2, data, L=NULL, B=10, weight_optimiser="Training MSE")
  expect_type(tst, "list")
})

test_that("crf works for: Pointwise variance; AR(1); L>20", {
  tst <- crf(y~x.1+x.2, data, L=40, B=10, weight_optimiser="Pointwise variance", correlation="ar1", x0=data.frame(x.1=1,x.2=1))
  expect_type(tst, "list")
})
