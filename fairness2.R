library(foreach)
library(doParallel)

##helper functions
#point_estimate
pe = function(x, p) {
  n = length(x)
  mu_h = mean(x)
  
  q_h = quantile(x, probs = p)
  m_h = mean(ifelse(x > q_h, 0, x)) / mu_h
  
  return(m_h)
}

#bootstrap
bts = function (x , b = 500, p) {
  res = rep (NA, b)
  n = length (x)
  for (ii in 1 : b) {
    xb = sample(x , size = n , replace = TRUE)
    qtb = quantile(xb, probs = p)
    yb = ifelse(xb > qtb, 0, xb)
    res[ii] = mean(yb) / mean(xb)
  }
  return (var(res))
}

#closed-form
analytical = function(x, p) {
  n = length(x)
  mu_h = mean(x)
  
  q_h = quantile(x, probs = p)
  m_h = mean(ifelse(x > q_h, 0, x)) / mu_h
  
  y = x * (ifelse(x <= q_h, 1, 0) - m_h)
  z = ifelse(x <= q_h, 1, 0) - p
  
  #first accounts for uncertainty of q_hat, second doesn't
  return(list(var_h_q_h = mean((y - q_h * z) ^ 2) / (n * mu_h ^ 2), 
              var_h_q_fixed = mean(y ^ 2) / (n * mu_h ^ 2)
              )
         )
}

#for simulations
simu_lnorm_once = function(n, p, meanlog, sdlog) {
  x = rlnorm(n = n, meanlog = meanlog, sdlog = sdlog)
  q_h = quantile(x, probs = p)
  
  #(closed-form) variance estimates
  var_h = analytical(x = x, p = p)
  return(list(m_h = mean(ifelse(x > q_h, 0, x)) / mean(x), 
              var_h_q_h = var_h$var_h_q_h,
              var_h_q_fixed = var_h$var_h_q_fixed, 
              var_h_boot = bts(x = x, p = p)))
}

#parallel computing
simu_lnorm = function(n_iterations = 5000, n, p = 0.75, meanlog, sdlog, n_clusters = 10) {
  #initiate cluster for parallel computing
  cluster = makeCluster(n_clusters)
  clusterExport(cluster, c('pe', 'bts', 'analytical', 'simu_lnorm_once'))
  registerDoParallel(cluster)
  
  res = foreach(ii = 1:n_iterations, .combine = rbind) %dopar% {
    tmp = simu_lnorm_once(n = n, p = p, meanlog = meanlog, sdlog = sdlog)
    data.frame(tmp)
  }
  
  #stop cluster
  stopCluster(cl = cluster)
  
  return(res)
}


##simulation studies
#cases
p = 0.75
mus = c(.4, -.3, .6)
sigmas = c(.5, 1, .5)
ns = c(2000, 5000, 10000)

res = data.frame(matrix(0, nrow = 9, ncol = 11))
colnames(res) = c('mu', 'sigma', 'n', 'm', 'var', 
                  'bias_h', 'bias_h_q_fix', 'bias_h_boot',
                  'cover_h', 'cover_h_q_fix', 'cover_h_boot')
res_raw = list()

counter = 0
for (ii in 1:3) {
  for (jj in 1:3) {
    counter = counter + 1
    print(paste('Now simulating the ', counter, '-th case!', sep = ''))
    
    #record data generation parameters
    res[counter, 1] = mus[ii]
    res[counter, 2] = sigmas[ii]
    res[counter, 3] = ns[jj]
    
    tmp = simu_lnorm(meanlog = mus[ii], sdlog = sigmas[ii], n = ns[jj])
    res_raw[[counter]] = tmp
    
    #compute ground truth
    q = qlnorm(p = p, meanlog = mus[ii], sdlog = sigmas[ii])
    #compute m from Butler and McDonald (1989)
    res[counter, 4] = plnorm(q = q, meanlog = mus[ii] + sigmas[ii] ^ 2, sdlog = sigmas[ii])
    #compute true variance
    res[counter, 5] = var(tmp$m_h)
      
    #compute relative biases
    res[counter, 6] = 100 * (mean(tmp$var_h_q_h) / res[counter, 5] - 1)
    res[counter, 7] = 100 * (mean(tmp$var_h_q_fixed) / res[counter, 5] - 1)
    res[counter, 8] = 100 * (mean(tmp$var_h_boot) / res[counter, 5] - 1)
      
    #compute coverage rates
    res[counter, 9] = 100 * mean(tmp$var_h_q_h >= (tmp$m_h - res[counter, 4]) ^ 2 / qnorm(.975) ^ 2)
    res[counter, 10] = 100 * mean(tmp$var_h_q_fixed >= (tmp$m_h - res[counter, 4]) ^ 2 / qnorm(.975) ^ 2)
    res[counter, 11] = 100 * mean(tmp$var_h_boot >= (tmp$m_h - res[counter, 4]) ^ 2 / qnorm(.975) ^ 2)
  }
}
print(res)


##empirical example
library(AER)
data('CPS1988')
urb = CPS1988$wage[CPS1988$smsa == 'yes']
sbb = CPS1988$wage[CPS1988$smsa == 'no']
p = .75

print(c(group_size = length(urb), m_hat = pe(x = urb, p = p), var = analytical(x = urb, p = p)))
print(c(group_size = length(sbb), m_hat = pe(x = sbb, p = p), var = analytical(x = sbb, p = p)))