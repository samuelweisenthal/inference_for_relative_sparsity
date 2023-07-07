est.v0 = function(beta,eps,gamma,u){
  n=length(eps)
  V=0
  T=length(eps[[1]]$Ss)
  for (i in 1:n){
    # r = \sum_t a[t]s[t][2,]
    # ER = E \sum_t s[t][2,]E[A[t]|s[t]]
    
    for (t in 1:T){
      V = V - u^(t-1)*expit(beta%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]][2,]
    }
 
    }
  V/n -gamma*sum(beta**2)
}

