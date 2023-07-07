# simulates data
# can simulate discrete and continuous a
# can also compute rewards if desired as an observation

library(mvtnorm)
source('utils.R')


gen.data = function(b0,n,init.state.mean,init.state.var,trans.var,T,tau,R){
  
  eps = list()
  for (i in 1:n){
    
    s = t(mvtnorm::rmvnorm(n=1,mean=init.state.mean,
                           sigma=init.state.var))
    
    As = list()
    Ss = list()
    Rs = list()
    
    for (t in 1:T){
      Ss[[t]] = s
      
      a = rbinom(1,1,prob=expit(t(b0)%*%s))
      s = t(mvtnorm::rmvnorm(n=1,mean=s+matrix(tau*rep(a,dim(s)[1])),
                             sigma=trans.var))
      
      As[[t]] = a
      Rs[[t]] = R(Ss[[t]],As[[t]],T)
    }
    eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
  
}

sim.db = function(b0,N,tau,T,R,sd.epsi,mean.first,sd.first,plot.=FALSE,seed=1){
  #sd.epsi = 2
  # obsY = matrix(nrow=N,ncol=T)
  # obsA = matrix(nrow=N,ncol=T)
  # obsMU = matrix(nrow=N,ncol=T+1)
  set.seed(seed)
  K = dim(b0)[1]
  eps = list()

  for (i in 1:N){
    Y = list()#rep(NA,T)
    A = list()#rep(NA,T)
    MU = list()#rep(NA,T+1)
    rr = list()#rep(NA,T)
    mm = list()
    MU[[1]] = matrix(rep(1,K))
    Y[[1]] = t(rmvnorm(1,mean.first,sd.first))
    A[[1]] = rbinom(1,1,expit(t(b0)%*%Y[[1]]))
    rr[[1]] = R(Y[[1]],A[[1]],T)
    mm[[1]] = 0
    if (T>1){
      #MU[[2]] = MU[[1]]*(1-tau*A[[1]]+tau*(1-A[[1]]))
      MU[[2]] = MU[[1]]*(1+tau*A[[1]]) # works now. w multivar? yes.
      for (t in 2:T){
        # epsi = rnorm(1,0,sd.epsi)
        epsi = t(rmvnorm(1,matrix(rep(0,K)),sd.epsi**2))
        #y =(Y[[t-1]] - MU[[t-1]] + epsi)/((1+sd.epsi**2)**(1/2)) + MU[[t]]
        Y[[t]] = sqrt(solve(sd.epsi**2+diag(K)))%*%(Y[[t-1]] - MU[[t-1]] + epsi) + MU[[t]]
        #Y[[t]] = matrix(y)
        A[[t]] = rbinom(1,1,expit(t(b0)%*%Y[[t]]))
        MU[[t+1]] = MU[[t]]*(1+tau*A[[t]])
        mm[[t]] = (Y[[t]]>0)*1
        #MU[[t+1]] = MU[[t]]*(1-tau*A[[t]]+tau*(1-A[[t]]))
        
        # uncomment to censor.
        censor=0
        if (censor){
        if (mm[[t]] == 1){
          Y[[t]] = NA
        }
        }
        rr[[t]] = R(Y[[t]],A[[t]],T)
        #rr[[t]] = R(Y[[t]],A[[t]],Y[[t+1]],T)
      }
    }
    # obsY[i,] = unlist(Y)
    # obsA[i,] = unlist(A)
    # obsMU[i,]=unlist(MU)
    eps[[i]] = list(As=A,Ss=Y,Rs=rr,Ms=mm)
    
  }
  if (plot.){
    par(mfrow=c(2,2))
    plot(1:T,rep(0,T),ylim=c(min(obsY),max(obsY)),type='l',lty=2)
    for (i in 1:N){
      lines(1:T,obsY[i,])
    }
    
    plot(1:T,apply(obsY,2,var),type='l')
    plot(1:T,apply(obsY,2,mean),type='l')
    plot(1:(T+1),apply(obsMU,2,mean),type='l')
  }
  eps
}


sim.db.nocen = function(b0,N,tau,T,R,sd.epsi,mean.first,sd.first,plot.=FALSE){
  #sd.epsi = 2
  # obsY = matrix(nrow=N,ncol=T)
  # obsA = matrix(nrow=N,ncol=T)
  # obsMU = matrix(nrow=N,ncol=T+1)
  K = dim(b0)[1]
  eps = list()
  
  for (i in 1:N){
    Y = list()#rep(NA,T)
    A = list()#rep(NA,T)
    MU = list()#rep(NA,T+1)
    rr = list()#rep(NA,T)
    mm = list()
    MU[[1]] = matrix(rep(1,K))
    Y[[1]] = t(rmvnorm(1,mean.first,sd.first))
    A[[1]] = rbinom(1,1,expit(t(b0)%*%Y[[1]]))
    rr[[1]] = R(Y[[1]],A[[1]],T)
    mm[[1]] = 0
    if (T>1){
      #MU[[2]] = MU[[1]]*(1-tau*A[[1]]+tau*(1-A[[1]]))
      MU[[2]] = MU[[1]]*(1+tau*A[[1]]) # works now. w multivar? yes.
      for (t in 2:T){
        # epsi = rnorm(1,0,sd.epsi)
        epsi = t(rmvnorm(1,matrix(rep(0,K)),sd.epsi**2))
        #y =(Y[[t-1]] - MU[[t-1]] + epsi)/((1+sd.epsi**2)**(1/2)) + MU[[t]]
        Y[[t]] = sqrt(solve(sd.epsi**2+diag(K)))%*%(Y[[t-1]] - MU[[t-1]] + epsi) + MU[[t]]
        #Y[[t]] = matrix(y)
        A[[t]] = rbinom(1,1,expit(t(b0)%*%Y[[t]]))
        MU[[t+1]] = MU[[t]]*(1+tau*A[[t]])
        mm[[t]] = (Y[[t]]>0)*1
        #MU[[t+1]] = MU[[t]]*(1-tau*A[[t]]+tau*(1-A[[t]]))
        
        # uncomment to censor.
        censor=0
        if (censor){
          if (mm[[t]] == 1){
            Y[[t]] = NA
          }
        }
        rr[[t]] = R(Y[[t]],A[[t]],T)
      }
    }
    # obsY[i,] = unlist(Y)
    # obsA[i,] = unlist(A)
    # obsMU[i,]=unlist(MU)
    eps[[i]] = list(As=A,Ss=Y,Rs=rr,Ms=mm)
    
  }
  if (plot.){
    par(mfrow=c(2,2))
    plot(1:T,rep(0,T),ylim=c(min(obsY),max(obsY)),type='l',lty=2)
    for (i in 1:N){
      lines(1:T,obsY[i,])
    }
    
    plot(1:T,apply(obsY,2,var),type='l')
    plot(1:T,apply(obsY,2,mean),type='l')
    plot(1:(T+1),apply(obsMU,2,mean),type='l')
  }
  eps
}

# variable length
sim.db.random.T = function(b0,N,tau,T,R,sd.epsi,mean.first,sd.first){
  
  K = dim(b0)[1]
  Ts = sample(1:T,n,replace=TRUE)
  eps = list()
  for (i in 1:N){
    Y = list()#rep(NA,T)
    A = list()#rep(NA,T)
    MU = list()#rep(NA,T+1)
    rr = list()#rep(NA,T)
    mm = list()
    MU[[1]]=matrix(rep(1,K))
    Y[[1]] = t(rmvnorm(1,mean.first,sd.first))
    A[[1]]=rbinom(1,1,expit(t(b0)%*%Y[[1]]))
    rr[[1]] = R(Y[[1]],A[[1]],T)
    mm[[1]] = 0
    if (Ts[i]>1){
      #MU[[2]] = MU[[1]]*(1-tau*A[[1]]+tau*(1-A[[1]]))
      MU[[2]] = MU[[1]]*(1+tau*A[[1]]) # works now. w multivar? yes.
      for (t in 2:Ts[i]){
        epsi = t(rmvnorm(1,matrix(rep(0,K)),sd.epsi**2))
        Y[[t]] = sqrt(solve(sd.epsi**2+diag(K)))%*%(Y[[t-1]] - MU[[t-1]] + epsi) + MU[[t]]
        A[[t]] = rbinom(1,1,expit(t(b0)%*%Y[[t]]))
        MU[[t+1]] = MU[[t]]*(1+tau*A[[t]])
        mm[[t]] = (Y[[t]]<0)*1
        rr[[t]] = R(Y[[t]],A[[t]],T)
        
      }
    }
    
    eps[[i]] = list(As=A,Ss=Y,Rs=rr,Ms=mm)
    
  }
  eps
}
# sim.episodes = function(nep,eplengths,S,A,actions,pib,P_s,P_sasp,
#                         phi,basis.degree,seed,par,intercept,init.state.var,
#                         init.state.par,trans.par,trans.var,disc.a,dim.s,R){
#   set.seed(seed)
#   eps = list()
#   for (i in 1:nep){
#     
#     s = t(mvtnorm::rmvnorm(n=1,mean=init.state.par,
#                            sigma=init.state.var))
#     #s=exp(s)
#     eplen = eplengths[i]
#     
#     As = list()
#     Ss = list()
#     Rs = list()
#     
#     for (t in 1:eplen){
#       Ss[[t]] = s
#       
#       if(disc.a){
#         dist.a = pib(s=s,a=A,par=par,
#                      phi=phi,basis.degree=basis.degree,
#                      intercept=intercept)
#         
#         a = sample(A,1,prob=dist.a)
#       }else{
#         pol.par = get.cont.policy.mean(s,a,par,phi,basis.degree,intercept)
#         pol.mean = pol.par$mean
#         pol.sd = pol.par$sd
#         
#         #if (fix.sim.sd.1){pol.sd=1}
#         
#         a = rnorm(n=1,mean=pol.mean,sd=pol.sd) 
#       }
#       
#       if (eplen>1){
#         
#         if (disc.a){
#           smean = get.transit.mean(a,s,trans.par,actions,dim.s)
#         }else{
#           smean = get.transit.mean.cont.a(a,s,trans.par,dim.s)
#         }
#         
#         s = t(mvtnorm::rmvnorm(n=1,mean=smean,
#                                sigma=trans.var))
#         #s=exp(s)
#       }
#       As[[t]] = a
#       Rs[[t]] = get.R.by.arg(R,Ss[[t]],As[[t]],s)
#     }
#     eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
#   }
#   eps
# }
