rm(list=ls())

typen = "III"
epsilon_u = 1

## if you want the so-called 'error function'
erf <- function(x) {2.0 * pnorm(x * sqrt(2.0)) - 1.0}
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'


fun_F_x <- function(x) {
  return(.5-.5 * erf((log(-x+1.0)+.5*log(2.0))/sqrt(2.0*log(2.0))))
}
# integrate of F(x)
fun_int_F<-function(a,b) {
  integrand <- function(x) {1.0-fun_F_x(x)}
  temp = integrate(integrand,lower = a, upper = b)
  return(temp$value)
}