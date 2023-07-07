source('is.R')

mi.ind = function(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u,sigmas,sel){
  #psa(beta,s,a)/psa(b,s,a)*r-gamma*beta**2
  beta=c(beta1,beta2)
  b=c(b1,b2)

  mi(beta,b,Ss,As,gamma,Rs,u,sigmas,sel)
}

#####


# outer product here bc s is Dx1, a-exp is 1x1.
dthetalogpsa = function(theta,eta,s,a,sel){
#dthetalogpsa = function(theta,eta,s,a,sel.p,sel.s){
  #note, this doesn't work for dbdlogpsa. need the sel diff
  if(is.na(matrix(s)[1,])){
    0
  }else{
    #(s*sel)%*%(a-ps1(theta,eta,s,sel))
    (a-ps1(theta,eta,s,sel))*(s*sel)
  }
} 

detalogpsa = function(theta,eta,s,a,sel){
  #new 
  #(s*(1-sel))%*%(a-ps1(theta,eta,s,sel))
  (a-ps1(theta,eta,s,sel))*(s*(1-sel))
  
}

detalogpsa.b = function(eta,s,a){
  #new 
  #(s*(1-sel))%*%(a-ps1(theta,eta,s,sel))
  (a-ps1.b(eta,s))*s
  
}

#dthetaexpit = function(x){expit(x)*(1-expit(x))}
dthetaexpit = function(theta,eta,s,a,sel){
  lp = t(theta)%*%(s*sel)+t(eta)%*%(s*(1-sel))
  expit(lp)*(1-expit(lp))
}

dthetaexpit.b =  function(eta,s,a){
  lp = t(eta)%*%s
  expit(lp)*(1-expit(lp))
}

ddthetalogpsa = function(theta,eta,s,a,sel){
  # theta, eta given together and picked out by sel
  if(is.na(matrix(s)[1,])){
    0
  }else{
  
    #-1*dthetaexpit(t(theta)%*%(s*sel)+t(eta)%*%(s*(1-sel)))*(s*sel)%*%t(s*sel)
    -1*dthetaexpit(theta,eta,s,a,sel)*(s*sel)%*%t(s*sel)
  }
}

detadthetalogpsa = function(theta,eta,s,a,sel){
  #new
    #-1*dthetaexpit(t(theta)%*%(s*sel)+t(eta)%*%(s*(1-sel)))*(s*sel)%*%t(s*(1-sel))
    -1*dthetaexpit(theta,eta,s,a,sel)*(s*sel)%*%t(s*(1-sel))
}

dthetalogkappa = function(theta,eta,Ss,As,sel){
  dlk=0
  T=length(As)
  for (t in 1:T){
    dlk=dlk+dthetalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  dlk
}

ddthetalogkappa = function(theta,eta,Ss,As,sel){
  dlk=0
  T=length(As)
  for (t in 1:T){
    dlk=dlk+ddthetalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  dlk
}

detadthetalogkappa = function(theta,eta,Ss,As,sel){
  # new
  # is it right?
  dedtlogk = 0 
  T=length(As)
  for (t in 1:T){
    dedtlogk = dedtlogk + detadthetalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  dedtlogk
}

dthetakappa = function(theta,eta,Ss,As,sel){
  kappa(theta,eta,Ss,As,sel)*dthetalogkappa(theta,eta,Ss,As,sel)
}


detakappa = function(theta,eta,Ss,As,sel){
  # not used
  #new 
  # correct???? yes, same as above...
  kappa(theta,eta,Ss,As,sel)*detalogkappa(theta,eta,Ss,As,sel)
}

detakappa.b = function(eta,Ss,As){
  # not used
  #new. it s d eta kappa_b,b. so no beta. no need for sel
  kappa.b(eta,Ss,As)*detalogkappa.b(eta,Ss,As)
}

detalogkappa = function(theta,eta,Ss,As,sel){
  # new
  # correct?
  delogk = 0
  T=length(As)
  for (t in 1:T){
    delogk = delogk + detalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  delogk
}

detalogkappa.b = function(eta,Ss,As){
  # new
  # correct?
  delogk = 0
  T=length(As)
  for (t in 1:T){
    delogk = delogk + detalogpsa.b(eta,Ss[[t]],As[[t]])
  }
  delogk
}

ddthetakappa = function(theta,eta,Ss,As,sel){
  #new
  (dthetakappa(theta,eta,Ss,As,sel)%*%t(dthetalogkappa(theta,eta,Ss,As,sel))
  + kappa(theta,eta,Ss,As,sel)*ddthetalogkappa(theta,eta,Ss,As,sel))
}
#ad.dthetakappa=Deriv(kappa,x='theta')


dbetari = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){
  
  dthetakappa(beta,b,Ss,As,sel)/kappa(b,b,Ss,As,sel)
}

dbetavi = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){

  G(Rs,u)*dbetari(beta,b,Ss,As,Rs,gamma,u,sigmas,sel)
}


dbri = function(beta,b,Ss,As,Rs,u,sel){
  
  (detakappa(beta,b,Ss,As,sel)/kappa.b(b,Ss,As) 
   + (-1)*kappa(beta,b,Ss,As,sel)/(kappa.b(b,Ss,As)^2)*detakappa.b(b,Ss,As))
}

dbvi = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){
  
  G(Rs,u)*dbri(beta,b,Ss,As,Rs,u,sel)
}


dbetapen=function(gamma,sigmas,beta,b,Ss,As,sel,other.args){
  # use.l sel here?

  if (other.args$use.l){
    K = dim(Ss[[1]])[1]
    #-gamma*dli(Ss,As,beta,rep(1,K))
    -gamma*dli.pen(Ss,As,beta,b,sel)
  }else{
  2*gamma*sigmas*(beta-b)
  }
}

dbetami = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel,other.args){
  # 2gamma is such a sledgehammer. need more specialized penalty
  
  # note that this only works for beta-b now...can 't change it..
  # should ideally have function for each piece, ie debtaV+gamma*dbetapen
  
  #G(Rs,u)*dthetakappa(beta,b,Ss,As,sel)/kappa(b,b,Ss,As,sel)-2*gamma*sigmas*(beta-b)
  
  dbetavi(beta,b,Ss,As,Rs,gamma,u,sigmas,sel) - dbetapen(gamma,sigmas,beta,b,Ss,As,sel,other.args)
}


ddbetari = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){
  (1/kappa(b,b,Ss,As,sel))*(dthetakappa(beta,b,Ss,As,sel)%*%t(dthetalogkappa(beta,b,Ss,As,sel)) 
                        + kappa(beta,b,Ss,As,sel)*ddthetalogkappa(beta,b,Ss,As,sel))
}

ddbetavi =  function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){
  G(Rs,u)*ddbetari(beta,b,Ss,As,Rs,gamma,u,sigmas,sel)
  
  }


ddbetapen = function(gamma,sigmas,beta,b,Ss,As,sel,other.args){
  # use.l sel here?

  if (other.args$use.l){
    K = dim(Ss[[1]])[1]
    #-gamma*ddli(Ss,As,beta,rep(1,K))
    -gamma*ddli.pen(Ss,As,beta,b,sel)
  }else{
  2*gamma*diag(length(sigmas))*sigmas
  }
}

ddbetami = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel,other.args){
  
  # I did algebra just using indices to check. you just put kappa into the vector
  # then you can use standard product rule on the entries of the vector.
  
  #G(Rs,u)/kappa(b,b,Ss,As,sel)*(dthetakappa(beta,b,Ss,As,sel)%*%t(dthetalogkappa(beta,b,Ss,As,sel)) 
  #                        + kappa(beta,b,Ss,As,sel)*ddthetalogkappa(beta,b,Ss,As,sel)) - 2*gamma*diag(length(sigmas))*sigmas

  ddbetavi(beta,b,Ss,As,Rs,gamma,u,sigmas,sel) - ddbetapen(gamma,sigmas,beta,b,Ss,As,sel,other.args)
  }


dbdbetari = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){
  #dthetakappa(beta,Ss,As)%*%t(dthetakappa(b,Ss,As) ) != dthetakappa(beta,Ss,As)%*%t(dthetakappa(b,Ss,As) )
  # since eg db1dbeta2 != db2 dbeta1
  
  K = dim(beta)[1]

  ((-1)*(kappa(b,b,Ss,As,sel)^(-2))*dthetakappa(beta,b,Ss,As,sel)%*%t(dthetakappa(b,b,Ss,As,rep(1,K))) #b^2
    + (1/kappa(b,b,Ss,As,sel))*dthetakappa(beta,b,Ss,As,sel)%*%t(detalogkappa(beta,b,Ss,As,sel))
    + kappa(beta,b,Ss,As,sel)/kappa(b,b,Ss,As,sel)*detadthetalogkappa(beta,b,Ss,As,sel))
}


dbdbetavi = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel){

  G(Rs,u)*dbdbetari(beta,b,Ss,As,Rs,gamma,u,sigmas,sel)
}

dbdbetapen = function(gamma,sigmas,beta,b,Ss,As,sel,other.args){
  # use.l sel here?
 
  if(other.args$use.l){
    #new
    K = dim(Ss[[1]])[1]
    #-dbdbetali(Ss,As,beta,b,rep(1,K))
    -dbdbetali(Ss,As,beta,b,sel)
  }else{
  (-1)*2*gamma*diag(length(sigmas))*sigmas
  }
}

dbdbetami = function(beta,b,Ss,As,Rs,gamma,u,sigmas,sel,other.args){

  K = dim(beta)[1]
  
  dbdbetavi(beta,b,Ss,As,Rs,gamma,u,sigmas,sel) - dbdbetapen(gamma,sigmas,beta,b,Ss,As,sel,other.args)
}


dli = function(Ss,As,theta,sel){
  # note that this just treats variable length episodes
  # as variable length.  I wonder actually, do we need to impute
  # if we have variable length episodes - eg with dropout.
  dl=0
  T=length(Ss)
  for (t in 1:T){
    dl=dl+dthetalogpsa(theta,theta,Ss[[t]],As[[t]],sel)
  }
  dl
}

dli.pen = function(Ss,As,theta,eta,sel){
  # note that this just treats variable length episodes
  # as variable length.  I wonder actually, do we need to impute
  # if we have variable length episodes - eg with dropout.
  dl=0
  T=length(Ss)
  for (t in 1:T){
    dl=dl+dthetalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  dl
}

dlit = function(s,a,theta){
  # where is this used?
    dthetalogpsa(theta,s,a)
	}

ddli = function(Ss,As,theta,sel){
  ddl=0 
  T=length(Ss)
  for (t in 1:T){
    ddl=ddl+ddthetalogpsa(theta,theta,Ss[[t]],As[[t]],sel)
    }
  ddl
}

ddli.pen = function(Ss,As,theta,eta,sel){
  ddl=0 
  T=length(Ss)
  for (t in 1:T){
    ddl=ddl+ddthetalogpsa(theta,eta,Ss[[t]],As[[t]],sel)
  }
  ddl
}

dbdbetali = function(Ss,As,theta,eta,sel){
  #new
  dbdbetal=0 
  T=length(Ss)
  for (t in 1:T){
    #detadthetalogpsa = function(theta,eta,s,a,sel)
    dbdbetal=dbdbetal+detadthetalogpsa(theta=theta,eta=eta,Ss[[t]],As[[t]],sel)
  }
  dbdbetal
}

dbdbetali.pen = function(Ss,As,theta,eta,sel){
  #new
  dbdbetal=0 
  T=length(Ss)
  for (t in 1:T){
    #detadthetalogpsa = function(theta,eta,s,a,sel)
    dbdbetal=dbdbetal+detadthetalogpsa(theta=theta,eta=eta,Ss[[t]],As[[t]],sel)
  }
  dbdbetal
}


dbdbetami.pd = function(beta,b,Ss,As,Rs,gamma,u){
  T = length(As)
  dbdbetav = 0
  for (t in 1:T){
    dbdbetav=dbdbetav + (-1)*u^(t-1)*Rs[[t]]*(kappa(b,Ss[1:t],As[1:t])^(-2))*dthetakappa(beta,Ss[1:t],As[1:t])%*%t(dthetakappa(b,Ss[1:t],As[1:t])) 
  }
  dbdbetav + 2*gamma*diag(length(b))
}

mi.ind.pd = function(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u){
  #psa(beta,s,a)/psa(b,s,a)*r-gamma*beta**2
  beta=c(beta1,beta2)
  b=c(b1,b2)
  
  mi.pd(beta,b,Ss,As,gamma,Rs,u)
  
}

dbetami.pd = function(beta,b,Ss,As,Rs,gamma,u){
  # 2gamma is such a sledgehammer. need more specialized penalty
  
  # note that this only works for beta-b now...can 't change it..
  # should ideally have function for each piece, ie debtaV+gamma*dbetapen
  T = length(As)
  dbetav = 0
  for (t in 1:T){
    ss = Ss[1:t]
    aa = As[1:t]
    dbetav = dbetav + dthetakappa(beta,ss,aa)/kappa(b,ss,aa)*u^{t-1}*Rs[[t]]
  }
  
  dbetav-2*gamma*(beta-b)
}

ddbetami.pd = function(beta,b,Ss,As,Rs,gamma,u){
  
  ddbetav = 0
  T=length(As)
  for (t in 1:T){
    
    ddbetav = ddbetav + u^(t-1)*Rs[[t]]/kappa(b,Ss[1:t],As[1:t])*(dthetakappa(beta,Ss[1:t],As[1:t])%*%t(dthetalogkappa(beta,Ss[1:t],As[1:t])) 
                                                                  + kappa(beta,Ss[1:t],As[1:t])*ddthetalogkappa(beta,Ss[1:t],As[1:t])) 
  }  
  ddbetav-2*gamma*diag(length(beta))
}

#mi = function(beta,b,s,a,gamma){
#  psa(beta,s,a)/psa(b,s,a)*R(s,a)#-gamma*beta**2
#}

# dbetam = function(beta,b,eps,gamma,u){
#   # where used?
#   n=length(eps)
#   db = 0
#   for (i in 1:n){
#     db = db+dbetami(beta,b,Ss,As,Rs,gamma,u)
#   }
#   db/n
# }


