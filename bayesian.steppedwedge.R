
## Bayesian Hierarchical inference by JAGS
bayes.pwr = function(xxx, fam, p, t=6, Nij=NULL, N=20, tx.eff, mu0, m=0, sigma.e=1,
                     sigma.a=0.5, alpha=0.05, tprec, sig.0=10,
                     sig.e.s=0.0001, sig.e.r=0.0001, sig.a.s=0.0001, sig.a.r=0.0001, model.file=NULL) {
  # fam = outcome type ("gaussian", "poisson", or "binomial")
  # p = number of clusters
  # t = no. of time periods (baseline period included)
  # Nij = matrix of # of at-risk individuals by cluster(i)-period(j)
  # N = equal # obs in each cluster-period, used if Nij is NULL
  # tx.eff = treatment effect (RD for gaussian, RR for poisson, OR for binomial)
  # mu0 = individual level baseline mean response (natural scale)
  # m = range of time trend coef as a fraction of the intervention effect coef
  # (+ve m corresponds to increasing coef over time)
  # sigma.e = individual level SD (for gaussian outcome)
  # sigma.a = SD of cluster means (linear predictor scale)
  # alpha = significance level
  # tprec = the time effect precision in JAGS script (1/variance)
  # sig.0 = variance of prior for a0 ("intercept" term, linear predictor scale)
  # sig.e.s = shape parameter of gamma hyperprior for sigma.e
  # sig.e.r = rate parameter of gamma hyperprior for sigma.e
  # sig.a.s = shape parameter of gamma hyperprior for random effects
  # sig.a.r = rate parameter of gamma hyperprior for random effects
  # model.file = (optional) name of file specifying analytic model
  fams = c("gaussian", "poisson", "binomial") # set allowed outcome types
  if (!(fam %in% fams)) return(
    "Allowed outcome distributions are \"gaussian\", \"poisson\", or \"binomial\".")
  if (is.null(Nij)) Nij = matrix(N, nrow=p, ncol=t)
  # if no model file is specified, create name of default file from outcome type
  if (is.null(model.file)) mf = paste(fam, "txt", sep=".")
  
  # number of iterations for inference ##
  n_iter = 5000
  # design matrix (cluster level)
  x = matrix(1, nrow=p, ncol=t)
  tt = (0:(p-1)) %% (t-1) + 1 # transition times
  for (i in 1:p) x[i, 1:tt[i]] = 0
  
  # cluster effects
  a = rnorm(p, 0, sigma.a)
  a0 = switch(fam, "gaussian"=mu0, "binomial"=log(mu0 / (1-mu0)),
              "poisson"=log(mu0))
  theta = switch(fam, "gaussian"=tx.eff, "binomial"=log(tx.eff),
                 Page 1
                 10.7_Supplementary_Code (2).txt
                 "poisson"=log(tx.eff))
  b = m * theta * (0:(t-1)) / (t-1)
  ### simulating data ###
  dat = matrix(NA, nrow=p, ncol=t)
  if (fam=="gaussian") {
    sigma.ep = sigma.e / sqrt(Nij)
    for (i in 1:p) {
      mu.i = a0 + a[i] + b + theta * x[i, ]
      dat[i, ] = rnorm(t, mu.i, sigma.ep[i, ])
    }
  }
  if (fam=="poisson") {
    for (i in 1:p) {
      lambda.i = exp(log(Nij[i, ]) + a0 + a[i] + b + theta * x[i, ])
      dat[i, ] = rpois(t, lambda.i)
    }
  }
  if (fam=="binomial") {
    for (i in 1:p) {
      pi.i = 1 / (1 + exp(-(a0 + a[i] + b + theta * x[i, ])))
      dat[i, ] = rbinom(t, Nij[i, ], pi.i)
    }
  }
  ### calling JAGS in R ###
  jags_works = 0; attempt = 1
  while (jags_works == 0){
    # Following call should work will with all outcome types; not all parameters passed
    # are needed for every outcome type
    jags=try(jags.model(file = mf,
                        data = list("t" = t, "p" = p, "y" = dat, "x" = x, "N" = Nij, "tprec" =
                                      tprec,
                                    "sig.0" = sig.0, "sig.e.s" = sig.e.s, "sig.e.r" = sig.e.r,
                                    "sig.a.s" = sig.a.s, "sig.a.r" = sig.a.r),
                        inits = list("a" = rep(0,p), "b" = rep(0,t), "theta.x" = 0, "a0" = 1,
                                     "sig.e" = 1, "prec" = 0.5),
                        quiet=T,n.chains = 2) )
    jags_works = as.numeric(length(jags) > 1 | attempt == 10)
    attempt = attempt + 1
  }
  update(jags, 1500) # more burn-in
  my.samples = coda.samples(jags, variable.names = c("a","b","theta.x"), n.iter =
                              n_iter)
  check.this = coda.samples(jags, variable.names = c("sigma.a"), n.iter = n_iter)
  samples = my.samples[[1]]
  samples2 = my.samples[[2]]
  l = as.numeric(quantile(samples[, (p+t+1)], c(alpha/2))) # 95% BIC
  Page 2
  10.7_Supplementary_Code (2).txt
  u = as.numeric(quantile(samples[, (p+t+1)], c(1-alpha/2)))
  R.h0 = ifelse(u < 0 | l > 0, 1, 0) # change one-sided to 2-sided
  if (fam=="binomial" | fam=="poisson"){
    trt.effect = exp(samples[, (p+t+1)]) # one chain
    trt.effect2 = exp(samples2[, (p+t+1)])
    thetax=mcmc(trt.effect)
    thetax2=mcmc(trt.effect2)
    comb=mcmc.list(thetax,thetax2)
    diag_MCMC=gelman.diag(comb) # Gelman-Rubin diagnostics
    r_hat=as.numeric(unlist(diag_MCMC)[1]) # the point estimate
    ess=effectiveSize(trt.effect) # effective sample size
    out = c(mean(trt.effect), R.h0, exp(l), exp(u),
            mean(check.this[[1]]),r_hat,ess,trt.effect=list(trt.effect),trt.effect2=list(trt.eff
                                                                                         ect2)) # check.this[[1]] is the variance of the random effect
    class(out)='trial'
  }
  if (fam=="gaussian"){
    trt.effect = samples[, (p+t+1)]
    trt.effect2 = samples2[, (p+t+1)]
    thetax=mcmc(trt.effect)
    thetax2=mcmc(trt.effect2)
    comb=mcmc.list(thetax,thetax2)
    diag_MCMC=gelman.diag(comb) # Gelman-Rubin diagnostics
    r_hat=as.numeric(unlist(diag_MCMC)[1]) # the point estimate
    ess=effectiveSize(trt.effect)
    out = c(mean(trt.effect), R.h0, l, u,
            mean(check.this[[1]]),r_hat,ess,trt.effect=list(trt.effect),trt.effect2=list(trt.eff
                                                                                         ect2))
    class(out)='trial'
  }
  
  return(out)
}

##### End of the main function #####
### JAGS model for count outcome - save as a file with name "poisson.txt" ###
model{
  for ( j in 1:p ){ # cluster
    for ( k in 1:t ){ # time period
      y[j,k] ~ dpois(lambda[j,k])
      log(lambda[j,k]) <- log(N[j,k]) + a0 + a[j] + b[k] + x[j,k] * theta.x
    } }
  # priors
  theta.x ~ dnorm(0, 1/9)
  a0 ~ dnorm(0, 1/sig.0^2) # baseline individual level outcome
  Page 3
  10.7_Supplementary_Code (2).txt
  for (i in 1:p){ # cluster effect
    a[i] ~ dnorm(0, prec) }
  for (i in 1:t){ # fixed effect
    b[i] ~ dnorm(0, tprec) }
  prec ~ dgamma(sig.a.s, sig.a.r)
  sigma.a <- 1 / prec
}
### JAGS model for binary outcome - save as a file with name "binomial.txt" ###
model{
  for ( j in 1:p ){ # cluster
    for ( k in 1:t ){ # time period
      y[j,k] ~ dbinom(lambda[j,k], N[j,k])
      logit(lambda[j,k]) <- a0 + a[j] + b[k] + x[j,k]*theta.x
    } }
  # priors
  theta.x ~ dnorm(0, 1/9)
  a0 ~ dnorm(0, 1/sig.0^2) # baseline log-odds for an individual outcome
  for (i in 1:p){ # cluster effect
    a[i] ~ dnorm(0 , prec) }
  for (i in 1:t){ # fixed effect
    b[i] ~ dnorm(0 , tprec) }
  prec ~ dgamma(sig.a.s, sig.a.r)
  sig.e ~ dgamma(sig.e.s, sig.e.r)
  sigma.a <- 1/prec
}
### JAGS model for continuous outcome - save as a file with name "gaussian.txt" ###
model{
  for ( j in 1:p ){ # cluster
    for ( k in 1:t ){ # time period
      y[j,k] ~ dnorm(mu[j,k], N[j,k] / sig.e^2)
      mu[j,k] <- a0 + a[j] + b[k]+ x[j,k] * theta.x
    } }
  # priors
  theta.x ~ dnorm(0, 1/9)
  a0 ~ dnorm(0, 1/sig.0^2) # baseline mean outcome
  for (i in 1:p){ # cluster effect
    a[i] ~ dnorm(0 , prec) }
  for (i in 1:t){ # fixed time effect
    b[i] ~ dnorm(0 , tprec) }
  prec ~ dgamma(sig.a.s, sig.a.r)
  sig.e ~ dgamma(sig.e.s, sig.e.r)
  sigma.a <- 1/prec
}
