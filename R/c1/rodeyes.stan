data{
int<lower=0> N;
int<lower=0> y;
real<lower=0> aprior;
real<lower=0> bprior;
}
parameters{
real<lower=0,upper=1> theta;
}
model{
y ~ binomial(N,theta);
theta ~ beta(aprior, bprior);
}
