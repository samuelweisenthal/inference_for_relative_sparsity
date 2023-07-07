source('grad.R')
source('is.R')
library(MASS)
psi=1
if (psi){
  solve=ginv
}
##
get.var = function(eps,beta,b,gamma,u,sigmas,dont.est.beh=FALSE,sel,weighted.var,other.args){
  
  n = length(eps)
  
  # it's weird, these are usually vectors or matrices, but I init them as 0. Must broadcast?
  DDM = 0
  DDL = 0
  cross = 0
  
  DDBETAVN = 0
  RN = 0 
  DBETAVN = 0
  DBETARN = 0
  DDBETARN = 0
  VN = 0
  
  DBDBETAVN = 0
  DBRN = 0
  DBVN = 0
  DBDBETARN = 0
  
  DDBETAPEN = 0
  DBDBETAPEN = 0
  
  for (i in 1:n){
    # if (is.na(eps[[i]]$Ss[[2]])){
    #   eps[[i]]$Ss[[2]] = 99
    # }
  
    DDM = DDM+ddbetami(beta=beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,gamma,u,sigmas,sel,other.args)
    cross = cross+dbdbetami(beta=beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,gamma,u,sigmas,sel,other.args)
    # we give DDL no selection, we just care about est of b
    DDL = DDL + ddli(eps[[i]]$Ss,eps[[i]]$As,theta=b,sel=rep(1,dim(beta)[1]))
    
    DDBETAVN = DDBETAVN + ddbetavi(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    RN = RN + ri(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,sel=sel)
    DBETAVN = DBETAVN + dbetavi(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    DBETARN = DBETARN + dbetari(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    DDBETARN = DDBETARN + ddbetari(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    VN = VN + vi(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,u=u,sel=sel)
    
    DBDBETAVN = DBDBETAVN +  dbdbetavi(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    DBRN = DBRN + dbri(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,u=u,sel=sel)
    DBVN = DBVN + dbvi(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    DBDBETARN = DBDBETARN + dbdbetari(beta=beta,b=b,Ss=eps[[i]]$Ss,As=eps[[i]]$As,Rs=eps[[i]]$Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
    
    # use.l sel here?
    DDBETAPEN = DDBETAPEN + ddbetapen(gamma,sigmas,beta,b,eps[[i]]$Ss,eps[[i]]$As,sel,other.args)
    DBDBETAPEN = DBDBETAPEN + dbdbetapen(gamma,sigmas,beta,b,eps[[i]]$Ss,eps[[i]]$As,sel,other.args)
  }
  

  DDM = DDM/n
  cross = cross/n
  DDL = DDL/n
  
  DDBETAVN =DDBETAVN/n
  RN = RN/n
  DBETAVN = DBETAVN/n
  DBETARN = DBETARN/n
  DDBETARN = DDBETARN/n
  VN = VN/n
  
  DBDBETAVN = DBDBETAVN/n
  DBRN = DBRN/n
  DBVN = DBVN/n
  DBDBETARN = DBDBETARN/n
  
  DDBETAPEN = DDBETAPEN/n
  DBDBETAPEN = DBDBETAPEN/n
  
  # -ddbetapen(gamma,sigmas,beta,b,Ss,As,sel) needs to be eg average!!! need to build above
  #ddbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  #dbdbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  # I think minus sign ok in front of penalty deriv? check. yes - all d*pen are negative, but we want them all positive

  # old, mistake
  #DDM.W = DDBETAVN/RN + (-1)*(DBETAVN/RN)%*%t(DBETARN/RN) + (-1)*(DBETARN/RN)%*%t(DBETAVN/RN) + (-1)*(VN/RN)*(DDBETARN/RN) + ((-2)*(-1)*(VN/RN)*(DBETARN/RN)*2*RN)%*%t(DBETARN/RN) - DDBETAPEN#ddbetapen(gamma,sigmas,beta,b,Ss,As,sel)

  DDM.W = DDBETAVN/RN + (-1)*(DBETAVN/RN)%*%t(DBETARN/RN) + (-1)*(DBETARN/RN)%*%t(DBETAVN/RN) + (-1)*(VN/RN)*(DDBETARN/RN) + ((-2)*(-1)*(VN/RN)*(DBETARN/RN))%*%t(DBETARN/RN) - DDBETAPEN#ddbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  
  
  #DDM.W = DDBETAVN/RN + (-1)*DBETAVN%*%t(DBETARN)/(RN^2) + (-1)*DBETARN%*%t(DBETAVN)/(RN^2) + (-1)*VN*DDBETARN/(RN^2)+((-2)*(-1)*(VN)*DBETARN*2*RN)%*%t(DBETARN)/(RN^3) - DDBETAPEN#ddbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  
  #print("weighted hessian")
  #print(DDM.W)
  #DDM.W = apply((DDM.W),c(1,2),function(x){min(x,1e30)})
  
  # old, mistake
  #cross.w = DBDBETAVN/RN + (-1)*(DBETAVN/RN)%*%t(DBRN/RN) + (-1)*(DBETARN/RN)%*%t(DBVN/RN) + (-1)*(VN/RN)*(DBDBETARN/RN) + ((-1)*(-2)*(VN/RN)*(DBETARN/RN)*2*RN)%*%t(DBRN/RN) - DBDBETAPEN#dbdbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  
  cross.w = DBDBETAVN/RN + (-1)*(DBETAVN/RN)%*%t(DBRN/RN) + (-1)*(DBETARN/RN)%*%t(DBVN/RN) + (-1)*(VN/RN)*(DBDBETARN/RN) + ((-1)*(-2)*(VN/RN)*(DBETARN/RN))%*%t(DBRN/RN) - DBDBETAPEN#dbdbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  #browser()
  
  #cross.w = DBDBETAVN/RN + (DBETAVN/(RN^2)*(-1))%*%t(DBRN) + (-1)*DBETARN%*%t(DBVN)/(RN^2) + (-1)*VN*DBDBETARN/(RN^2) +((-1)*(-2)*VN*DBETARN*2*RN)%*%t(DBRN)/(RN^3) - DBDBETAPEN#dbdbetapen(gamma,sigmas,beta,b,Ss,As,sel)
  

  VAR = matrix(0,nrow=length(beta),ncol=length(beta))
  VAR.B = matrix(0,nrow=length(b),ncol=length(b))
  
  VAR.W = matrix(0,nrow=length(beta),ncol=length(beta))
  MEAN.W = matrix(0,nrow=length(beta),ncol=1)
  
  if (dont.est.beh){
    # because these terms are multiplied by sqrtn(bn-b0), which is zero if bn=b0
    cross = 0
    cross.W = 0
  }
  
  for(i in 1:n){
    dbetam = dbetami(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,gamma,u,sigmas,sel,other.args)
    #dl = dli(eps[[i]]$Ss,eps[[i]]$As,b,sel)
    dl.b = dli(eps[[i]]$Ss,eps[[i]]$As,b,sel=rep(1,dim(beta)[1])) # actaully, isn't this one unecessary?
    ct = solve(DDL)%*%dl.b
    # KxK KX1 - KxK KxK Kx1, where K = dim(state) 
    
    H = -solve(DDM)%*%(dbetam - cross%*%ct)
    H.B = -solve(DDL)%*%dl.b
    VAR = VAR + H%*%t(H)
    VAR.B = VAR.B + H.B%*%t(H.B)
    
    # is the penalty term right. do i need to divide by n or something or take sum?
    # penalty deriative is in P_n. And - is correct, because we want it positive
    # use.l sel here?
    dbetam.w = dbetavi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,gamma,u,sigmas,sel)/RN + (-1)*vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u,sel)*DBETARN/(RN^2)-dbetapen(gamma,sigmas,beta,b,eps[[i]]$Ss,eps[[i]]$As,sel,other.args)
    H.W = -solve(DDM.W)%*%(dbetam.w - cross.w%*%ct)
    VAR.W = VAR.W + H.W%*%t(H.W)
    MEAN.W = MEAN.W + H.W
    
  } 
  
  
  VAR = VAR/n
  VAR.B = VAR.B/n
  
  MEAN.W = MEAN.W/n
  MEAN.W.sq = MEAN.W%*%t(MEAN.W)

  VAR.W = VAR.W/n
  
  sub.mean.sq=1
  if (sub.mean.sq){
    VAR.W = VAR.W - MEAN.W.sq
  }

  if (weighted.var){
    # note that will still return unweighted DDM,cross. These aren't ever used
    # in calling script.
    VAR=VAR.W
  }
  
  list(DDM=DDM,cross=cross,DDL=DDL,VAR=VAR,VAR.B=VAR.B)
}

# should be something like
emp.exp = function(eps,beta,b,gamma,u,f){

n=length(eps)
Ts = unlist(lapply(eps,function(x){length(x$As)}))
Tcal = unique(Ts)

el=0
for (tp in Tcal){
  # which episodes have length tp
  Itp=which(Ts==tp)
  # consider that set of episodes
  epstp = eps[Itp]
  ntp = length(epstp)
  
  Seltp = 0
  for (ep in epstp){
    
    Seltp = Seltp+f(beta=beta,b,ep$Ss,ep$As,ep$Rs,gamma,u)
    
  }
  
  Eeltp = (1/ntp)*Seltp
 
  Eel = Eel + (ntp/n)*Eeltp
}

Eel

}

get.var.randomT = function(eps,beta,b,gamma,u,dont.est.beh=FALSE){

  n=length(eps)
  Ts = unlist(lapply(eps,function(x){length(x$As)}))
  Tcal = unique(Ts)
  
  EDDM = 0
  EDDL = 0
  Ecross = 0
  for (tp in Tcal){
    # which episodes have length tp
    Itp=which(Ts==tp)
    # consider that set of episodes
    epstp = eps[Itp]
    ntp = length(epstp)
    
    SDDMtp = 0
    SDDLtp = 0
    Scrosstp = 0
    for (ep in epstp){
      
      SDDMtp = SDDMtp+ddbetami(beta=beta,b,ep$Ss,ep$As,ep$Rs,gamma,u)
      Scrosstp = Scrosstp+dbdbetami(beta=beta,b,ep$Ss,ep$As,ep$Rs,gamma,u)
      #print(cross)
      SDDLtp = SDDLtp + ddli(ep$Ss,ep$As,theta=b)
      
    }
    
    EDDMtp = (1/ntp)*SDDMtp
    Ecrosstp = (1/ntp)*Scrosstp
    EDDLtp = (1/ntp)*SDDLtp
    
    EDDM = EDDM + (ntp/n)*EDDMtp
    EDDL = EDDL + (ntp/n)*EDDLtp
    Ecross = Ecross + (ntp/n)*Ecrosstp
  }
  
  
  EVAR=matrix(0,nrow=length(beta),ncol=length(beta))
  EVAR.B=matrix(0,nrow=length(b),ncol=length(b))
  
  for (tp in Tcal){
    # which episodes have length tp
    Itp=which(Ts==tp)
    # consider that set of episodes
    epstp = eps[Itp]
    ntp=length(epstp)
    
    SVARtp = 0
    SVAR.Btp = 0
    for (ep in epstp){
      
      dbetam=dbetami(beta,b,ep$Ss,ep$As,ep$Rs,gamma,u)
      dl = dli(ep$Ss,ep$As,b)
      # DxD DX1 - DxD DxD Dx1 
      
      ct = -solve(EDDL)%*%dl
      H = -solve(EDDM)%*%(dbetam + Ecross%*%ct)
      H.B = -solve(EDDL)%*%dl
      
      SVARtp = SVARtp + H%*%t(H)
      SVAR.Btp = SVAR.Btp + H.B%*%t(H.B)
      
    }
    
    EVARtp = (1/ntp)*SVARtp
    EVAR.Btp = (1/ntp)*SVAR.Btp
    
    EVAR = EVAR + (ntp/n)*EVARtp
    EVAR.B = EVAR.B + (ntp/n)*EVAR.Btp
  }
  
 
  list(DDM=EDDM,cross=Ecross,DDL=EDDL,VAR=EVAR,VAR.B=EVAR.B)
}