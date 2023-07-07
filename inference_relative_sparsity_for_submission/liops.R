# computes terms in log likelihood that depend on policy

# ok I think to comment, not used anywhere
#logpsa = function(theta,eta,s,a,sel){log(psa(theta,eta,s,a,sel))}

li = function(b,Ss,As,sel){
  l=0;
  T=length(Ss);
  for (t in 1:T){
    l=l+slog(psa(b,b,Ss[[t]],As[[t]],sel))
    };
  l
  }

l = function(b,eps,sel){
  ll=0
  for (ep in eps){
    ll=ll+li(b,ep$Ss,ep$As,sel)
  }
  ll
}


li.pen = function(beta,b,Ss,As,sel){
  l=0;
  T=length(Ss);
  for (t in 1:T){
    l=l+slog(psa(beta,b,Ss[[t]],As[[t]],sel))
  };
  l
}

l.pen = function(beta,b,eps,sel){

  K = dim(eps[[1]]$Ss[[1]])[1] 
  ll=0
  for (ep in eps){
    ll=ll+li.pen(beta,b,ep$Ss,ep$As,sel) - li.pen(b,b,ep$Ss,ep$As,rep(1,K))
  }
  ll
}

l.nt = function(b,eps){
  # compute likelihood treating all nT as independent
  # not used.  just to test against other way (n obs)
  l=0
  T=length(eps[[1]]$As)
  for (ep in eps){
    for (t in 1:T){
      l = l + li(b,ep$Ss[t],ep$As[t])
    }
  }
  l
}
#####
# pre-new scripts

l.mdp = function(P_s,P_sasp,pi_sa,Ss,As,S,A,par,phi,basis.degree,intercept,
                 init.state.par,init.state.var,tilt,conf.control,use.nonpara.conf.contr){
  
  # NOTE: this isn't the real likelihood bc I ignore terms 
  # that dnd on the policy, bc the max of policy dnd on them

  pi_sa.val = pi_sa(s=Ss[[1]],a=As[[1]],par=par,phi=phi,
                    basis.degree=basis.degree,intercept=intercept)
  # likelihood dnd
  # uncomment p_s here and below in l expression to get true likelihood
  #p_s.val = P_s(Ss[[1]],init.state.par,init.state.var)#get.p(Ss[[1]],p=P_s)
  l = eps.log(pi_sa.val) #+ eps.log(p_s.val) 
  
  epL =  length(Ss)
  if(epL>1){
    # if one stage, just use thing above
  for(i in 1:(epL-1)){ # normally need first. I guess this is ok for check

    
    ##
    # was
    #l = l + eps.log(P_sasp[Ss[i],As[i],Ss[i+1]]) + eps.log(pi_sa[Ss[i],As[i]])
    ##
    pi_sa.val = pi_sa(s=Ss[[i+1]],a=As[[i+1]],par=par,phi=phi,
                      basis.degree=basis.degree,intercept=intercept)
    # uncomment below and in l to get true likelihood
    #p_sasp.val = P_sasp(Ss[[i]],As[[i]],Ss[[i+1]])
    l = l + eps.log(pi_sa.val) #+ eps.log(p_sasp.val) 
  }
  }
  l
}


ltotal.mdp = function(P_s,P_sasp,pi_sa,episodes,S,A,
                      par,phi,basis.degree,intercept,
                      gamma=NULL,R=NULL,pib_sa=NULL,parb=NULL,tilt=FALSE,conf.control,
                      use.nonpara.conf.contr){

  # computes likelihood of set of mdp episodes
  ltot = 0
  for (episode in episodes){
    ltot = ltot + l.mdp(P_s=P_s,P_sasp=P_sasp,pi_sa=pi_sa,
                        Ss=episode$Ss,As=episode$As,S=S,A=A,
                        par=par,phi=phi,basis.degree=basis.degree,
                        intercept=intercept,tilt=tilt,conf.control=conf.control,
                        use.nonpara.conf.contr=use.nonpara.conf.contr)
  }
  ltot
}
