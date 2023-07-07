# functions to center and scale data in time-dependent manner
# and put into different data structures
# code is slow, but it is not bottleneck

get.rep.stat = function(ar,stat){
  mean.ar = t(apply(X=ar,MARGIN=c(2,3),FUN=stat))
  ar.mean.rep = array(NA,dim=dim(ar))
  
  N = dim(ar)[1]
  T = dim(ar)[2]
  K = dim(ar)[3]
  for (k in 1:K){
    for (i in 1:N){
      #print(mean.ar[i,])
      ar.mean.rep[i,,k] = mean.ar[k,]
    }
  }
  ar.mean.rep
}

ep.to.arr = function(e){
  N=length(e)
  Ts = lapply(e,function(x){length(x$As)})
  T=max(unlist(Ts))
  K=dim(e[[1]]$Ss[[1]])[1]
  sa = array(NA,c(N,T,K))

  for (i in 1:N){
    for (j in 1:T){
      for (k in 1:K){

        sa[i,j,k] = e[[i]]$Ss[[j]][k,]
      }
    }
  }
  sa
}

arr.to.ep = function(arr,eps){
  N = length(eps)
  Ts = lapply(eps,function(x){length(x$As)})
  T=max(unlist(Ts))
  K=dim(eps[[1]]$Ss[[1]])[1]
  
  for (i in 1:N){
    for (j in 1:T){
      for (k in 1:K){
        eps[[i]]$Ss[[j]][k,]=arr[i,j,k]
      }
    }
  }
  eps
}

fill.in.states = function(e,sm.cs,dim.s){
  # put matrix of states into episode list object
  eplengths = unlist(lapply(e,function(x){length(x$As)}))
  k = 1
  for (j in 1:length(e)){
    for (i in 1:eplengths[j]){
      e[[j]]$Ss[[i]] = matrix(sm.cs[k,],nrow=dim.s,ncol=1)
      k=k+1
    }
  }
  e
}

scale.arr = function(eps){
  a=ep.to.arr(eps)
  av=get.rep.stat(a,var)
  a.s=a/sqrt(av)
  eps.s=arr.to.ep(a.s,eps)
  list(a.s=a.s,eps.s=eps.s)
}

arr.to.ep.HAL = function(N,T,K,arr,e){
  # merge with above
  for (i in 1:N){
    for (j in 1:T){
      e[[i]]$Ss[[j]]=matrix(NA,nrow=K,ncol=1)
      for (k in 1:K){
        e[[i]]$Ss[[j]][k,]=arr[i,j,k]
      }
    }
  }
  e
}

arr.to.mat = function(N,T,K,sa.cs){
  sa.cs.m = matrix(nrow=N*T,ncol=K)
  k=1
  for (i in 1:N){
    for (j in 1:T){
      sa.cs.m[k,]=sa.cs[i,j,]
      k=k+1
    }
  }
  sa.cs.m
}

center.scale.s = function(e,means=NULL,sds=NULL){
 # takes data as in sim.episodes()

  e.raw = e
  e.ti = e
  ss=lapply(e,function(x){x$Ss})
  dim.s=dim(e[[1]]$Ss[[1]])[1]
  eplengths = unlist(lapply(ss,function(x){length(x)}))
  #nsteps = sum(eplengths)
  # does it make sense to treat all time steps the same like this?
  # rather than eg scaling individually time step 1, time step 2, ...
  
  # moved to sm below
  #sm = matrix(unlist(ss),nrow=nsteps,ncol=dim.s,byrow=TRUE)

  ####################################################
  # test case
  run.test.case=0
  if (run.test.case){
  ar = array(seq(1,22,4),c(2,2,2))
  ar = array(runif(12),c(3,2,2))
  ar.mean.rep = get.rep.stat(ar,mean)
  #why doesn't sd work here. get sd=1. where did I set that?
  ar.sd.rep = get.rep.stat(ar,function(x){sqrt(var(x,na.rm=TRUE))})
  
  ar.c = ar - ar.mean.rep
  ar.cs = ar.c/ar.sd.rep
  }
  ####################################################
  
  N = length(e)
  T = eplengths[1]
  K = dim.s


  
  #if (T==1){
   
   
  sm = unlist.reshape.eps.ss.sumTsxK(ss)    
  #}else{
    #print("Need to debug scaler for Tneq1 for mimic")
   
#    sa = ep.to.arr(e) # why doesn't this work for mimic? 
 #   sm = arr.to.mat(N,T,K,sa)  #cleaner than above
    
  #    }


  if(is.null(means)){
    means = apply(sm,2,mean)
  }else{
    means=means
  }
  if(is.null(sds)){
    sds = apply(sm,2,function(x){sd(x,na.rm=TRUE)}) + 1e-20 
  }else{
    sds = sds
  }

  # note, no longer centering
  sm.c = sm#-tile(means,dim(sm))#t(replicate(dim(sm)[1],means))
  sm.cs = sm.c/tile(sds,dim(sm.c))#t(replicate(dim(sm)[1],sds))

  # fills in episode list object with scaled states
  e=fill.in.states(e,sm.cs,dim.s)

  list(e=e,means=means,sds=sds,
       sm.cs=sm.cs,sm.c=sm.c,sm=sm,
       e.raw=e.raw)
}

# center.scale.s = function(e,means=NULL,sds=NULL){
#   
#   e.raw = e
#   e.ti = e
#   ss=lapply(e,function(x){x$Ss})
#   dim.s=dim(e[[1]]$Ss[[1]])[1]
#   eplengths = unlist(lapply(ss,function(x){length(x)}))
#   #nsteps = sum(eplengths)
#   # does it make sense to treat all time steps the same like this?
#   # rather than eg scaling individually time step 1, time step 2, ...
#   
#   # moved to sm below
#   #sm = matrix(unlist(ss),nrow=nsteps,ncol=dim.s,byrow=TRUE)
#   
#   ####################################################
#   # test case
#   run.test.case=0
#   if (run.test.case){
#     ar = array(seq(1,22,4),c(2,2,2))
#     ar = array(runif(12),c(3,2,2))
#     ar.mean.rep = get.rep.stat(ar,mean)
#     #why doesn't sd work here. get sd=1. where did I set that?
#     ar.sd.rep = get.rep.stat(ar,function(x){sqrt(var(x))})
#     
#     ar.c = ar - ar.mean.rep
#     ar.cs = ar.c/ar.sd.rep
#   }
#   ####################################################
#   
#   N = length(e)
#   T = eplengths[1]
#   K = dim.s
#   
#   
#   sa = ep.to.arr(N,T,K,e)
#   sm = arr.to.mat(N,T,K,sa)  #cleaner than above
#   
#   
#   if(is.null(means)){
#     # if this is null then sd null too
#     sa.mean.rep = get.rep.stat(sa,mean)
#     means = apply(sm,2,mean)
#     # ok to divide by sd, not var?
#     # check for scale
#     sa.sd.rep = get.rep.stat(sa,function(x){sqrt(var(x))})
#     # latter in case sd =0. Never will happen but breaks code with small numbers of cases
#     # when using HAL I guess, since SK might be zero of one of the basis expansions.
#     # in general standardizing HAL not easy
#     # maybe no need to do so
#     sds = apply(sm,2,sd) + 1e-20 
#   }else{
#     sds = sds
#   }
#   
#   # note, no longer centering
#   sm.c = sm#-tile(means,dim(sm))#t(replicate(dim(sm)[1],means))
#   sm.cs = sm.c/tile(sds,dim(sm.c))#t(replicate(dim(sm)[1],sds))
#   sa.c = sa#-sa.mean.rep
#   sa.cs = sa.c/sa.sd.rep
#   
#   k = 1
#   for (j in 1:length(e)){
#     for (i in 1:eplengths[j]){
#       e[[j]]$Ss[[i]] = matrix(sm.cs[k,],nrow=dim.s,ncol=1)
#       k=k+1
#     }
#   }
#   
#   e.ti = arr.to.ep(N,T,K,sa.cs,e.ti)
#   
#   sa.cs.m = arr.to.mat(N,T,K,sa.cs)
#   # time sensitive is e.ti, sm.cs.ti
#   # non time sensitivie is e, sm.cs
#   sdotjk = 0
#   if (sdotjk){
#     e = e.ti
#     sm.cs = sa.cs.m
#   }
#   #plot(e[[1]],type='l')
#   #plot(sm[,2],type='l')
#   #plot(unlist(lapply())
#   list(e=e,means=means,sds=sds,
#        sm.cs=sm.cs,sm.c=sm.c,sm=sm,
#        e.raw=e.raw)
# }


#eps = sim.episodes(nep=100,eplengths=rep(2,100),S=S,A=A,actions=actions,
#                   pib=pi_sa,P_s=P_s,P_sasp=P_sasp,phi=phi,
#                   option=option,seed=1,par=par,
#                   intercept=intercept,init.state.var=init.state.var,
#                   init.state.par=init.state.par,trans.par=trans.par,
#                   trans.var=trans.var,
#                   disc.a=disc.a,dim.s=2)
#center.scale.s(eps,2)
