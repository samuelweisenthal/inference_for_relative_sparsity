#source('utils.R')
source('liops.R')
# expit defined in utils
expit = function(x){
    # k goes away, but stabilizes
    k  = max(0,x)#(0>x)*0+(x>=0)*x#max(0,x)
    as.vector(exp(x-k)/(exp(0-k)+exp(x-k)))
}

pen=function(beta,b,betahat,bhat,delta,sigmas){
  #print("Need to deal with adaptive lasso penalty")
  # maybe do need to scale in what
  #zou.gamma=1
  sum(what((betahat-bhat),delta)*abs(sigmas*(beta-b)))#^(zou.gamma)
}

gamm.pen = function(beta,b,sigmas){
  sum((sigmas*(beta-b))**2)
  #sum(beta**2)
}

G=function(Rs,u){
  g=0
  T=length(Rs)
  for (t in 1:T){
    g=g+u^(t-1)*Rs[[t]]
  }
  g
}

ps1 = function(theta,eta,s,sel){
  # issue with categorical cov. maybe after scaling?
  expit(t(theta)%*%(s*sel)+t(eta)%*%(s*(1-sel)))
}

ps1.b = function(eta,s){
  # issue with categorical cov. maybe after scaling?
  expit(t(eta)%*%s)
}

psa = function(theta,eta,s,a,sel){
  if(is.na(matrix(s)[1,])){
    1 # for censoring. not used now.
    }else{
  ps1(theta,eta,s,sel)^a*(1-ps1(theta,eta,s,sel))^(1-a)
  }
}

psa.b = function(eta,s,a){

    ps1.b(eta,s)^a*(1-ps1.b(eta,s))^(1-a)
}

kappa = function(theta,eta,Ss,As,sel){
  T=length(Ss)
  k=0
  for(t in 1:T){
    k=k+slog(psa(theta,eta,Ss[[t]],As[[t]],sel))
    }
  exp(k)
}

kappa.b = function(eta,Ss,As){
  T=length(Ss)
  k=0
  for(t in 1:T){
    k=k+slog(psa.b(eta,Ss[[t]],As[[t]]))
  }
  exp(k)
}


logkappa = function(theta,eta,Ss,As){
  slog(kappa(theta,eta,Ss,As))
}

logkappa.b = function(eta,Ss,As){
  slog(kappa.b(theta,eta,Ss,As))
}

# vi = function(beta,b,Ss,As,Rs,u){
#   kappab=1
#   kappabeta=1
#   T=length(Ss)
#   for(t in 1:T){
#     kappab=kappab*psa(b,Ss[[t]],As[[t]])
#     kappabeta = kappabeta*psa(beta,Ss[[t]],As[[t]])
#   }
#   G=G(Rs,u)
#   kappabeta/kappab*G  
# }


vi = function(beta,b,Ss,As,Rs,u,sel){
  log.kappab=0
  log.kappabeta=0
  T=length(As)
 # should use kappa function
  for(t in 1:T){
    log.kappab = log.kappab+slog(psa.b(b,Ss[[t]],As[[t]]))
    log.kappabeta = log.kappabeta+slog(psa(beta,b,Ss[[t]],As[[t]],sel))
  }
  G=G(Rs,u)
  
  exp(log.kappabeta)/exp(log.kappab)*G
}

vi.bfirst = function(b,beta,Ss,As,Rs,u,sel){
 vi(beta,b,Ss,As,Rs,u,sel)
}

gi = function(Rs,u){
  G(Rs,u)
}

ri = function(beta,b,Ss,As,sel){
  log.kappab=0
  log.kappabeta=0
  T=length(As)
  # should use kappa function
  for(t in 1:T){
    log.kappab = log.kappab+slog(psa.b(b,Ss[[t]],As[[t]]))
    log.kappabeta = log.kappabeta+slog(psa(beta,b,Ss[[t]],As[[t]],sel))
  }
  exp(log.kappabeta)/exp(log.kappab)
}

vi.spec = function(beta,b,Ss,As,Rs,u){
  kappab=1
  kappabeta=1
  T=length(As)
  suggs=rep(NA,T)
  behs=rep(NA,T)
  for(t in 1:T){
    beh=psa(b,Ss[[t]],As[[t]])
    kappab=kappab*beh
    behs[t]=beh
    sugg=psa(beta,Ss[[t]],As[[t]])
    kappabeta = kappabeta*sugg
    suggs[t]=sugg
  }
  G=G(Rs,u)
  list(v=kappabeta/kappab*G,rat=kappabeta/kappab,ret=G,
       kappabeta=kappabeta,kappab=kappab,suggs=round(suggs,2),behs=behs)
}


vn.randomT = function(beta,b,eps,u){
  # slow because function calls?
  # would be nice to make vi arg and then handle v0 this way, 
  # but too complicated
  n=length(eps)
  Ts = unlist(lapply(eps,function(x){length(x$As)}))
  Tcal = unique(Ts)
  
  EG = 0
  for (tp in Tcal){
    # which episodes have length tp
    Itp=which(Ts==tp)
    # consider that set of episodes
    epstp = eps[Itp]
    ntp=length(epstp)
    
    SGtp = 0
    for (ep in epstp){
      
      SGtp = SGtp +  vi(beta,b,ep$Ss,ep$As,ep$Rs,u)
      
    }
    
    EGtp = (1/ntp)*SGtp

    EG = EG + (ntp/n)*EGtp
  }
  EG
}

#UNWEIGHTED
weighted=1

if (weighted){
  vn = function(beta,b,eps,u,sel){
    n=length(eps)
    v=0
    w=0
    #ws = rep(NA,n)
    for (i in 1:n){
      v = v + vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u,sel)
      w = w + ri(beta,b,eps[[i]]$Ss,eps[[i]]$As,sel)
      #print(c("n=",n,"w=",w))
      #ws[i] = w
    }
    #print(c("v/n=",v/n,"v/w=",v/w))
    v/w
    #list(wvn=v/w,ws=ws)
  }
  
}else{

vn = function(beta,b,eps,u,sel){
 n=length(eps)
 v=0
 for (i in 1:n){
   v = v + vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u,sel)
 }
 v/n
}
}

# 



rn = function(beta,b,eps,u,sel){
  n=length(eps)
  ws = rep(NA,n)
  for (i in 1:n){
    ws[i] =  ri(beta,b,eps[[i]]$Ss,eps[[i]]$As,sel)
  }
  ws
}

sumrn = function(beta,b,eps,u,sel){
  sum(rn(beta,b,eps,u,sel))
}

#UNWEIGHTED
sigma.vn = function(beta,b,eps,u,sel){
 # try w mn
 vvn = vn(beta,b,eps,u,sel)
 ssumrn = sumrn(beta,b,eps,u,sel)
 n=length(eps)
 sigma2 = 0
 for (i in 1:n){
   sigma2 = sigma2 + (vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u,sel) - vvn)^2
 }
 # (sigma2/n)
 sqrt((1/(ssumrn/n)^2)*(sigma2/n))
}



# 
# sigma.vn = function(beta,b,eps,u,sel){
#   vvn = vn(beta,b,eps,u,sel)
#   ssumrn = sumrn(beta,b,eps,u,sel)
#   n=length(eps)
#   sigma2 = 0
#   for (i in 1:n){
#     sigma2 = sigma2 + (ri(beta,b,eps[[i]]$Ss,eps[[i]]$As,sel)/ssumrn)^2*(gi(eps[[i]]$Rs,u) - vvn)^2
#   }
#   # why no division by n?
#   sqrt(sigma2)
# }

sigma.mn = function(beta,b,eps,gamma,u,sigmas,sel){
  # try w mn
  mmn = mn(beta,b,eps,gamma,u,sigmas,sel)
  n=length(eps)
  sigma2 = 0
  for (i in 1:n){
    sigma2 = sigma2 + (mi(beta,b,eps[[i]]$Ss,eps[[i]]$As,gamma,eps[[i]]$Rs,u,sigmas,sel) - mmn)^2
  }
  
  sqrt(sigma2/n)
}


mi = function(beta,b,Ss,As,gamma,Rs,u,sigmas,sel,other.args){
  #psa(beta,s,a)/psa(b,s,a)*r-gamma*beta**2
 v.i=vi(beta,b,Ss,As,Rs,u,sel)

 K=dim(Ss[[1]])[1]

 if (other.args$use.l){
  # use.l sel here?  
 #pen=-gamma*li(beta,Ss,As,sel=rep(1,K))   
 pen=-gamma*li.pen(beta,b,Ss,As,sel=sel)  
 }else{
 pen=gamma*gamm.pen(beta,b,sigmas)
 }
 v.i-pen
}


# sums over t first, as would be derived typically.
# vn.pd = function(beta,b,eps,u){
#   T=length(eps[[1]]$As)
#   n=length(eps)
#   vn=0
#   for (t in 1:T){
#   vnt = 0
#   for (i in 1:n){
#     
#     log.kappab=0
#     log.kappabeta=0
#     for(k in 1:t){
#       log.kappab =       log.kappab+slog(psa(b,   eps[[i]]$Ss[[k]],eps[[i]]$As[[k]]))
#       log.kappabeta = log.kappabeta+slog(psa(beta,eps[[i]]$Ss[[k]],eps[[i]]$As[[k]]))
#     }
#     kappabeta=exp(log.kappabeta)
#     kappab=exp(log.kappab)
#     
#     vnt = vnt + kappabeta/kappab*eps[[i]]$Rs[[t]]*u^(t-1)
#       }
#   vn = vn + (1/n)*vnt
#   }
#   
# 
#   vn
# }

vi.pd = function(beta,b,Ss,As,Rs,u){
  vni = 0
  for (t in 1:T){
    
    log.kappab=0
    log.kappabeta=0
    for(k in 1:t){
      log.kappab =       log.kappab+slog(psa(b,   Ss[[k]],As[[k]]))
      log.kappabeta =    log.kappabeta+slog(psa(beta,Ss[[k]],As[[k]]))
    }
    kappabeta=exp(log.kappabeta)
    kappab=exp(log.kappab)
    
    vni = vni + kappabeta/kappab*Rs[[t]]*u^(t-1)
  }
  vni
}

# mi.pd = function(beta,b,Ss,As,gamma,Rs,u){
#   #psa(beta,s,a)/psa(b,s,a)*r-gamma*beta**2
#   v.i=vi.pd(beta,b,Ss,As,Rs,u)
#   v.i-gamma*gamm.pen(beta,b)
# }


# vn.pd = function(beta,b,eps,u){
#   T=length(eps[[1]]$As)
#   n=length(eps)
#   vn=0
#   for (i in 1:n){
#     vn = vn + vi.pd(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u)
#   }
#   vn/n
# }

# mn.pd = function(beta,b,eps,gamma,u,sigmas,pen.gam=TRUE){
#   if (pen.gam){
#     p = -gamma*gamm.pen(beta,b,sigmas)
#   }else{
#     p=l(beta,eps)
#   }
#   vn.pd(beta,b,eps,u) + p
# }

mn = function(beta,b,eps,gamma,u,sigmas,sel,other.args){
  use.l =other.args$use.l
  K=dim(eps[[1]]$Ss[[1]])[1]
  n = length(eps)
  
  if (!use.l){
    p = -gamma*gamm.pen(beta,b,sigmas)
  }else{
    # use.l sel here?
    #p=gamma*l(beta,eps,sel=rep(1,K))/n
    p=gamma*l.pen(beta,b,eps,sel=sel)/n
  }

   vn(beta,b,eps,u,sel) + p
}


# lambdan = function(n,delta){
#   #(2*log(n))^((1+delta)/2) #pg 1422
#   #log(n)
#   n^(1/100) # 
#   #0
# }

lambdan = function(n,delta){
  n^(1/100)
}
#lambdan = function(n,delta){n^(1-delta)}

what = function(D,delta){
  eps=1e-40
  1/(abs(D+eps)^(delta))
}

wn = function(beta,b,eps,gamma,u,betahat,bhat,delta,lambda,sigmas,sel,other.args){
  # maybe use LARS
  n = length(eps) #can rem
  mn(beta,b,eps,gamma,u,sigmas,sel,other.args) - lambda*pen(beta,b,betahat,bhat,delta,sigmas)
}

wn.pd = function(beta,b,eps,gamma,u,betahat,bhat,delta,lambda){
  # maybe use LARS
  n = length(eps) #can rem
  mn.pd(beta,b,eps,gamma,u) - lambda*pen(beta,b,betahat,bhat,delta)
}

v0i = function(beta,b,Ss,As,Rs,u,sel){
  # same for pd , should be
  T=length(As)
  V=0
  K=dim(Ss[[1]])[1]
  for (t in 1:T){
    if(K>1){
      if (is.na(matrix(Ss[[t]])[1,])){
      rt = 0
      }else{
      rt = -(expit(t(beta)%*%(Ss[[t]]*sel)+t(b)%*%(Ss[[t]]*(1-sel)))*Ss[[t]][2,])#/T
      }
    }else{
      if (is.na(matrix(Ss[[t]])[1,])){
        rt = 0
      }else{
        rt = -expit(t(beta)%*%(Ss[[t]]*sel)+t(b)%*%(Ss[[t]]*(1-sel)))*Ss[[t]]#/T        
      }
    }
    V = V + u^(t-1)*rt
  }
  V
  
}

est.m0.lambda = function(beta,eps,gamma,u,b,sigmas,sel,other.args){

  n=length(eps)
  V=0
  #K= dim(eps[[1]]$Ss[[1]])
  K=dim(eps[[1]]$Ss[[1]])[1]
  #T=length(eps[[1]]$Ss)
  for (i in 1:n){
    # r = \sum_t a[t]s[t][2,]
    # ER = E \sum_t s[t][2,]E[A[t]|s[t]]
    V = V + v0i(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u,sel)
  }
  

  if (!other.args$use.l){
    p = -gamma*gamm.pen(beta,b,sigmas)

  }else{
    # use.l sel here?
    #p = gamma*l(beta,eps,sel=rep(1,K))/n
    p = gamma*l.pen(beta,b,eps,sel=sel)/n
  }

  V/n + p
}

# est.m0.lambda.nocalls = function(beta,eps,gamma,u,b,pen.gam=TRUE){
#   # the pd and ordinary version of this are exactly equivalent
#   n=length(eps)
#   V=0
# 
#   for (i in 1:n){
#     
#     T=length(eps[[i]]$As)
#     Vi=0
#     #K=dim(eps[[i]]$Ss[[1]])[1]
#     
#     for (t in 1:T){
#       #if(K>1){
#         #if (is.na(eps[[i]]$Ss[[t]])){
#         #  rt = 0
#         #}else{
#        #   rt = -(expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]][2,])/T
#         #}
#     #  }else{
#         #if (is.na(eps[[i]]$Ss[[t]])){
#         #  rt = 0
#         #}else{
#           rt = -expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]]#/T        
#         #}
#      # }
#       Vi = Vi + u^(t-1)*rt
#     }
#     
#     V = V + Vi
#     
#     
#   }
#   if (pen.gam){
#     p = -gamma*gamm.pen(beta,b)
#   }else{
#     p=l(beta,eps)
#   }
#   
#   V/n + p
# }

# est.m0.lambda.nocalls.pd = function(beta,eps,gamma,u,b,pen.gam=TRUE){
#   # this is identical to the one above, the order of the sums are just changed.
#   n=length(eps)
#   T=length(eps[[1]]$As)
#   K=dim(eps[[1]]$Ss[[1]])[1]
#   
#   V=0  
#   for (t in 1:T){
#     
#     Vt=0
#     
#     for (i in 1:n){
#       if(K>1){
#         if (is.na(matrix(eps[[i]]$Ss[[t]])[1,])){
#           rt = 0
#         }else{
#           rt = -(expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]][2,])#/T
#         }
#       }else{
#         if (is.na(matrix(eps[[i]]$Ss[[t]])[1,])){
#           rt = 0
#         }else{
#           rt = -expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]]#/T        
#         }
#       }
#       Vt = Vt + u^(t-1)*rt
#     }
#     
#     V = V + Vt/n
#     
#     
#   }
#   if (pen.gam){
#     p = -gamma*gamm.pen(beta,b)
#   }else{
#     p=l(beta,eps)
#   }
#   
#   V + p
# }


#est.m0.lambda.nocalls(1,eps,0,u,1)

# est.m0.lambda.nocalls = function(beta,eps,gamma,u,b,pen.gam=TRUE){
#   
#   n=length(eps)
#   V=0
#   
#   for (i in 1:n){
#     
#     
#     T=length(eps[[i]]$As)
#     Vi=0
#     K=dim(eps[[i]]$Ss[[1]])[1]
#     
#     for (t in 1:T){
#       if(K>1){
#         if (is.na(eps[[i]]$Ss[[t]])){
#           rt = 0
#         }else{
#           rt = -(expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]][2,])/T
#         }
#       }else{
#         if (is.na(Ss[[t]])){
#           rt = 0
#         }else{
#           rt = -expit(t(beta)%*%eps[[i]]$Ss[[t]])*eps[[i]]$Ss[[t]]/T        
#         }
#       }
#       Vi = Vi + u^(t-1)*rt
#     }
#     
#     V = V + Vi
#     
#     
#   }
#   if (pen.gam){
#     p = -gamma*gamm.pen(beta,b)
#   }else{
#     p=l(beta,eps)
#   }
#   
#   V/n + p
# }

# est.m0.lambda.randomT = function(beta,eps,gamma,u,b,pen.gam=TRUE){
#   # note the same structure of mn. hence, could reduce code a lot
#   # by just making vi, v0i an argument.
#   n=length(eps)
#   Ts = unlist(lapply(eps,function(x){length(x$As)}))
#   Tcal = unique(Ts)
#   EG = 0
#   for (tp in Tcal){
#     # which episodes have length tp
#     Itp=which(Ts==tp)
#     # consider that set of episodes
#     epstp = eps[Itp]
#     ntp=length(epstp)
#     
#     SGtp = 0
#     for (ep in epstp){
#       SGtp = SGtp +  v0i(beta,ep$Ss,ep$As,ep$Rs,u)
#     }
#     
#     EGtp = (1/ntp)*SGtp
#     
#     EG = EG + (ntp/n)*EGtp
#   }
#   EG
#   
#   if (pen.gam){
#     p = -gamma*gamm.pen(beta,b)
#   }else{
#     p=l(beta,eps)
#   }
#   
#   EG + p
# }

# functions for (standard) importance sampling
# as a guide, used Thomas, Philip S. Safe reinforcement learning. Diss. 
# University of Massachusetts Libraries, 2015.


# ope.Thomas = function(episodes,pi_sa,par,parb,gamma,R,
#                       phi,basis.degree,intercept,
#                       S=NULL,A=NULL,P_s=NULL,P_sasp=NULL,
#                       conf.control,use.nonpara.conf.contr){
#   # specializing to binary action
#   # in the commented version above, it is more general, but slower
#   # uncomment all the functions above
#   
#   expit = function(x){
#     # k goes away, but stabilizes
#     k  = max(0,x)#(0>x)*0+(x>=0)*x#max(0,x)
#     exp(x-k)/(exp(0-k)+exp(x-k))
#   }
#   
#   beta=par
#   b=parb
#   n=length(episodes)
#   v=0
#   T=length(episodes[[1]]$Ss)
#   for (i in 1:n){
#     G=0
#     logkappa.beta=0
#     logkappa.b=0
#     for (t in 1:T){
#       st=episodes[[i]]$Ss[[t]]
#       rt=episodes[[i]]$Rs[[t]]
#       at=episodes[[i]]$As[[t]]
#       
#       G = G + gamma^(t-1)*rt
#       pi1.beta = expit(t(beta)%*%st)
#       logkappa.beta = logkappa.beta + log(pi1.beta^{at}*(1-pi1.beta)^(1-at))
#       
#       if (use.nonpara.conf.contr){
#         pi1.b=predict(conf.control,new_data=t(st))
#       }else{
#         pi1.b = expit(t(b)%*%st) 
#       }
#       logkappa.b = logkappa.b + log(pi1.b^{at}*(1-pi1.b)^(1-at))
#     }
#     v=v+exp(logkappa.beta)/exp(logkappa.b)*G
#   }
#   v=v/n
#   list(v=v,eps.info=NULL,var.v=NULL)
# }

# old=FALSE
# if (old){
# get.log.IS.ratio = function(ep,pi_sa,par,parb,t,phi,basis.degree,intercept,
#                             conf.control,use.nonpara.conf.contr){
# 
#   pie_sa.val = pi_sa(ep$Ss[[t]],ep$As[[t]],par,phi,basis.degree,intercept)
#   #print(ep$Ss[[t]])
#   pibe = pi_sa(ep$Ss[[t]],ep$As[[t]],parb,phi,basis.degree,intercept)
#   if (use.nonpara.conf.contr){
#     cc=predict(conf.control,new_data=t(ep$Ss[[t]]))
#     a=ep$As[[t]]
#     ccc = cc*(a==1)+(1-cc)*(a==0)
#     pib_sa.val =ccc
#   }else{
#     pib_sa.val = pibe
#   }
#   #print(c(pibe,ccc))
#   list(lisr=eps.log(pie_sa.val)-eps.log(pib_sa.val),pie_sa.val=pie_sa.val,
#        pib_sa.val=pib_sa.val)
# }
# 
# 
# rhoh.Thomas = function(ep,pi_sa,par,parb,gamma,R,phi,basis.degree,intercept,
#                        conf.control,use.nonpara.conf.contr){
# 
#   n.arg.R = length(formals(R))
#   ep.length = length(ep$As) #size of trajectory/episode
# 
#   if (n.arg.R==3){
#     ep.length = ep.length - 1
#   }
#   if (!is.function(R)){
#     R = ep$Rs
#   }
#   logm = 0
#   G = 0
#   ep.info = list()
#   for (t in 1:ep.length){
#     iso = get.log.IS.ratio(ep,pi_sa,par,parb,t,phi,basis.degree,intercept,
#                            conf.control,use.nonpara.conf.contr)
#     logisr = iso$lisr
#     Rt = get.Rt(ep,R,t)
#     logm = logm + logisr
#     # check against calcG. esp exp of gamma
#     G = G + (gamma^(t-1))*Rt
#     ep.info[[t]] = list(pie_sa.val=iso$pie_sa.val,pib_sa.val=iso$pib_sa.val)
#   }
#   ratio = exp(logm)
#   list(timed.info=ep.info,G=G,ratio=ratio,ratG = ratio*G)
#   #ratio*G
# }
# 
# ope.Thomas = function(episodes,pi_sa,par,parb,gamma,R,
#                       phi,basis.degree,intercept,
#                       S=NULL,A=NULL,P_s=NULL,P_sasp=NULL,
#                       conf.control,use.nonpara.conf.contr){
# 
#   IS = 0
# 
#   nD = length(episodes) # number episodes
#   eps.info = list()
#   j=1
#   for (ep in episodes){
# 
#     is.in =  rhoh.Thomas(ep,pi_sa,par,parb,
#                          gamma,R,phi,basis.degree,intercept,
#                          conf.control,use.nonpara.conf.contr)
#     IS = IS + is.in$ratG
#     j=j+1
#     eps.info[[j]] = is.in
#   }
#   v= IS/nD
#   var.v = 0
#   calc.var=0
#   if (calc.var){
#   for (ep in episodes){
#     is.in = rhoh.Thomas(ep,pi_sa,par,parb,
#                         gamma,R,phi,basis.degree,intercept,conf.control,use.nonpara.conf.contr)
#     var.v = var.v + (is.in$ratG-v)**2
#   }
#   # we divide once by n for approx the E
#   var.v = var.v/nD
#   # we divide again by n bc var v = (n*var v_i)/n^2
#   var.v = var.v/nD
#   }
#   list(v=v,eps.info=eps.info,var.v=var.v)
# }
# 
# 
# }
# 
# 
# 
# 
# 
# 
# 
