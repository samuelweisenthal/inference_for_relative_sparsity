library(latex2exp)
source('get.var.R')
source('plot.utils.R')
source('utils.R')

heightchunk=740
widthchunk = 2150
myheights=c(10,7,7)
myw = 21

my.cex.main=2.5 #main title
my.cex.axis=1.4 #numbers on axes
my.cex.lab=2.5 #labels on axes
cex.legend=1.3

#res.multiplier=100*1.5 # multiplied by ngam then
res.multiplier=100*1.5 # multiplied by ngam then

MC = function(M,n,b0,gamma,u,delta,init.state.mean,
              init.state.var,trans.var,T,minvar=FALSE,
              seed=1,tau,R,plt.obj,scale.s,sel,ii,other.args){
  print("need to make sim.db take other arguments")
  #set.seed(seed)
  exps = list()
  for (j in 1:M){
    if (j%%20 == 0) {
      print(c("MC",j))
    }
    
    #eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau,R)
    eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,sd.epsi=sd.epsi,mean.first=mean.first,
               sd.first=sd.first)
    lambda=lambdan(n,delta)
    exps[[j]]=process.eps(eps=eps,b0=b0,gamma=gamma,
                          u=u,delta=delta,lambda=lambda,plt.obj=plt.obj,
                          scale.s=scale.s,sel=sel,other.args=other.args)
    print(exps[[j]])
  }
  
  exps
}

process.eps = function(eps,b0,gamma,u,delta,lambda,plt.obj,scale.s,sel,split,other.args){
  #note b0 not used
  # scale.s 1 if you scale using sd. set sigmas for scale factors=1. 
  # otherwise if scale.s =0, use scale factors, set sigmas to sigmas from center.scale
  # too big. split up into pieces that process a dataset
  # and pieces necessary for MC. Maybe put MC pieces into MC function
  #so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K))
  
  n = length(eps)
  print(c("train n",n))
  K = dim(eps[[1]]$Ss[[1]])[1]
  T = length(eps[[1]]$Ss)
  
  hf = as.integer(n/2)
  n.test = NA#length(eps.test)
  #eps = eps[1:hf]
  n.train = NA#length(eps)
  
  if(split){
    # need to deal with this
    eps.test = eps[(hf+1):n]
    n.test = length(eps.test)
    print(c("n train sub test",n.test))
    
    # should call it eps.train
    eps = eps[1:hf]
    n.train = length(eps)
    print(c("n train sub train",n.train))
    
  }else{
    # just test and train on same set. when do I not split?
    # remember this is for calculating value on held out, not sample split
    # which is done elsewhere
    eps.test=eps
  }
  
  if(scale.s){
    #somehow this still seems buggy. check again
    so = center.scale.s(eps)
    sigmas = rep(1,K) # we don t do scale factor in penalty bc we scale state var
    eps.test = center.scale.s(eps.test,sds=so$sds)$e
    
  }else{
    # don t center and scale
    so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K)) #not centering and scaling
    # but we still need the sds for the scaling factors.
    soc = center.scale.s(eps)
    sigmas = soc$sds
    # generally, we would do the above with scale factors. instead
    # override scale factors by setting sigmas to 1. we are not scaling for simulations. 
    # we just generate data with variance 1
    #we just scale for real data
    sigmas = rep(1,K)
  }
  # a not affected by scaling
  aa = unlist(lapply(eps,function(x){x$As}))
  sm = so$sm.cs

  eps=so$e

  #eps.test=eps
  # how does 1:hf work if aa is dim nT
#  if (split){
    #bm=glm(aa[1:hf]~sm[1:hf,]-1,family='binomial')
    # because we made aa, sm on eps, which is epstrain. 
    # should have renamed eps to eps train sam
    # mistakes are never in the complex parts of the code
    # only in the parts that should be easy and you do quickly
    
 #   bm=glm(aa~sm-1,family='binomial') 
 # }else{
    bm=glm(aa~sm-1,family='binomial') 
#  }
  
  #a=d$a
  #s=d$s

  bn=bm$coefficients
  bn.pilot= bm$coefficients
  # to debug

  # set it again, since it is also used in variance calc
  #dont.est.beh=0 
  
  #Nelder-Mead gives convergence 10 degeneracy sometimes with K=5
  #if (K>1){
  #method='Nelder-Mead' # maybe write gradient descent using l1 approx
  #}else{
  method='BFGS' #must use BFGS for one dimensional case
  #}
  #browser()
  maxit=1e6
  est.bn.w.optim = 0
  #dont.est.beh=0
  if (other.args$dont.est.b){
    bn=b0
  }else{
  if (est.bn.w.optim){  
  
  # is it sel=rep(1,K)
  b.optim = optim(par=bn.pilot,l,
                 control=list(fnscale=-1,maxit=maxit),
                 method=method,
                 eps=eps,sel=rep(1,K))

  bn = b.optim$par
  }
  }

  #print(c("bn",bn))
  #start at behavioral to minimize positivity viol at beginning
  #print("commented out pd")
  #mn = mn.pd
  #vn = vn.pd

  #mn(rep(1,K),b=bn,eps=eps,gamma=gamma,u=u,sigmas=sigmas)
  
  
  optim.res = optim(par=c(bn),mn,
                    control=list(fnscale=-1,maxit=maxit),
                    method=method,
                    eps=eps,b=bn,gamma=gamma,u=u,sigmas=sigmas,sel=sel,other.args=other.args)
  

  #bn.optim = b.optim$par
  bn.optim=NA

  #bn.nt.optim = b2.optim$par
  #print(bn.nt.optim==bn.optim)
  #todo: check if optim w NT obs matches optim w N length T obs (check theory)
  #b=bn.optim
  betan=matrix(optim.res$par)
  # note it s on training sample
  rollout.betangamma=expit(sm%*%betan)
  rollout.bn = expit(sm%*%bn)
  pol.diff = abs(rollout.betangamma-rollout.bn)
  roll.diff.betangamma.bn.abs = sum(pol.diff)/(n*T)
  z = other.args$z

  roll.diff.p.gt.z  = sum(pol.diff>z)/(n*T)
  print(c("maxPolDiff",max(pol.diff)))
  vn.betagamma = vn(betan,bn,eps,u,sel)
  
  optim.res.lamb = optim(par=c(bn),wn,
                         control=list(fnscale=-1,maxit=maxit),
                         method=method,
                         b=bn,eps=eps,gamma=gamma,u=u,betahat=betan,
                         bhat=bn,delta=delta,lambda=lambda,sigmas=sigmas,sel=sel,other.args=other.args)
  print(c("converge?",optim.res.lamb$convergence))
  betanlambda=matrix(optim.res.lamb$par)
  rollout.betangammalambda = expit(sm%*%betanlambda)
  roll.diff.betangammalambda.bn.ab = sum(abs(rollout.betangammalambda-rollout.bn))/(n*T)
  # could rep vn with mn here?
 

  vn.betagammalambda=vn(betanlambda,bn,eps,u,sel)
  vn.betagammalambda.test = vn(betanlambda,bn,eps.test,u,sel)
  
  
  vn.bn=vn(bn,bn,eps,u,sel)
  var.v.bn=sigma.vn(bn,bn,eps,u,sel)
  
  #print(vn.betagammalambda>= vn.betagammalambda.test)
  wn.betagammalambda = wn(betanlambda,bn,eps,gamma,u,betan,bn,delta,lambda,sigmas,sel,other.args=other.args)
  wn.betagammalambda.test = wn(betanlambda,bn,eps.test,gamma,u,betan,bn,delta,lambda,sigmas,sel,other.args=other.args)
  
  tilde.optim.res = optim(par=c(bn),est.m0.lambda,
                          control=list(fnscale=-1,maxit=maxit),
                          method=method,
                          eps=eps,gamma=gamma,u=u,b=b0,sigmas=sigmas,sel=sel,other.args=other.args)
  
  tilde.betan=matrix(tilde.optim.res$par)
 
  var = get.var(eps,betan,bn,gamma,u=u,dont.est.beh=other.args$dont.est.b,sigmas=sigmas,sel=sel,weighted=other.args$weighted,other.args=other.args)$VAR
  var.plugintrue = get.var(eps,tilde.betan,b0,gamma,u=u,dont.est.beh=other.args$dont.est.b,sigmas=sigmas,sel=sel,weighted=other.args$weighted,other.args=other.args)$VAR
  var.betanlambda = get.var(eps,betanlambda,bn,gamma,u=u,dont.est.beh=other.args$dont.est.b,sigmas=sigmas,sel=sel,weighted=other.args$weighted,other.args=other.args)$VAR

  var.b = get.var(eps,betanlambda,bn,gamma,u=u,dont.est.beh=other.args$dont.est.b,sigmas=sigmas,sel=sel,weighted=other.args$weighted,other.args=other.args)$VAR.B
  
  sigma.vn. = sigma.vn(betanlambda,bn,eps,u,sel)
  sigma.vn.test = sigma.vn(betanlambda,bn,eps.test,u,sel)

  diff = betanlambda-bn
  
  
  # note use diff set in run.mc.R
  # hardcoded here for debugging, but not used otherwise
  use.diff = other.args$use.diff
  lamsel = abs(diff)>use.diff
  nsel=sum(lamsel)
  lp.rollout.betanlambda.sel = sm%*%(betanlambda*lamsel)
  lp.rollout.betanlambda.nosel = sm%*%(betanlambda*(1-lamsel))
  mean.lp.rollout.betanlambda.sel = mean(lp.rollout.betanlambda.sel)
  mean.lp.rollout.betanlambda.nosel = mean(lp.rollout.betanlambda.nosel)
  rollout.lps = expit(lp.rollout.betanlambda.sel+lp.rollout.betanlambda.nosel)[1:10]
  rollout.betangammalambda[1:10]
  comp.lp.sel.nosel = cbind(lp.rollout.betanlambda.sel,
                            lp.rollout.betanlambda.nosel)
  comp.lp.sel.nosel = comp.lp.sel.nosel[1:min(50,dim(comp.lp.sel.nosel)[1]),]
  
  ratios=rn(betanlambda,bn,eps,u,sel=rep(1,K))
  #print(c("betan",betan))
  #print("betanlambda")
  #print(betanlambda)
  #print("varbetanlambda")
  #print(sqrt(var.betanlambda/n.train)*1.96)
  list(var=var,
       var.betanlambda=var.betanlambda,
       var.b=var.b,
       var.plugintrue=var.plugintrue,
       gamma=gamma,
       delta=delta,
       lambda=lambda,
       lamsel=lamsel,
       diff.bn.v.betangammalambda=diff,
       lamsel=lamsel*1,
       nsel=nsel,
       lp.rollout.betanlambda.sel=lp.rollout.betanlambda.sel[1:100],
       mean.lp.rollout.betanlambda.sel=mean.lp.rollout.betanlambda.sel,
       lp.rollout.betanlambda.nosel=lp.rollout.betanlambda.nosel[1:100],
       mean.lp.rollout.betanlambda.nosel=mean.lp.rollout.betanlambda.nosel,
       comp.lp.sel.nosel=comp.lp.sel.nosel,
       betan=betan,
       betanlambda=betanlambda,
       bn=bn,
       bn.optim=bn.optim,
       tilde.betan=tilde.betan,
       vn.betagamma=vn.betagamma,
       vn.betagammalambda=vn.betagammalambda,
       vn.betagammalambda.test=vn.betagammalambda.test,
       vn.bn=vn.bn,
       var.v.bn=var.v.bn,
       wn.betagammalambda=wn.betagammalambda,
       wn.betagammalambda.test=wn.betagammalambda.test,
       sigam.vn. = sigma.vn.,
       sigam.vn.test = sigma.vn.test,
       split=split,
       n.train=n.train,
       n.test=n.test,
       rollout.betangamma =rollout.betangamma[1:100],
       rollout.betangammalambda = rollout.betangammalambda[1:100],
       rollout.betangamma.mean = mean(rollout.betangamma),
       rollout.betangammalambda.mean = mean(rollout.betangammalambda),
       rollout.bnmean = mean(rollout.bn),
       roll.diff.betangammalambda.bn.ab=roll.diff.betangammalambda.bn.ab,
       roll.diff.betangamma.bn.abs=roll.diff.betangamma.bn.abs,
       pol.diff=pol.diff[1:100],
       roll.diff.p.gt.z=roll.diff.p.gt.z,
       diff.bn.v.betangammalambda=diff,
       ratios=ratios)
}

lambda.path = function(eps.mimic,b0,gammas,lambdas,names,deltas,
                       file,sel,split,scale.s,lambda.select.ix=FALSE,
                       gamma.sel=1,other.args){
  
  # lambda path is misnomer, it iterates over delta, gamma, and lambda
  
  # should in some ways move split out of this function
  
  n=length(eps.mimic) # repeated!!! need to clean this up. get rid of inits in dim.R
  K = dim(eps.mimic[[1]]$Ss[[1]])[1]
  T = length(eps.mimic[[1]]$As)
  ti=Sys.time()
  outer.res=list()
  for (d in 1:length(deltas)){
    print(c("delta",deltas[d]))
    dres = list()
  for (j in 1:length(gammas)){
    print(c("gamma",gammas[j]))
    res = list()
    
    #nlam=2
    #par(mfrow=c(1,3))
    
    for (i in 1:length(lambdas)){
      print(c("lam",lambdas[i]))
      if ((gammas[j]>0) & (i==1)){plt.obj=1}else{plt.obj=0}
  
      res[[i]]=process.eps(eps.mimic,b0=b0,gamma=gammas[j],u=1,delta=deltas[d],
                           lambda=lambdas[i],plt.obj=0,scale.s=scale.s,sel=sel,
                           split=split,other.args=other.args)
    }

    dres[[j]] = res
  
  }
    outer.res[[d]] = dres
    
  }
  te = Sys.time()-ti
  print(te)
  # outer res is outer.res[d(delta)][j(gamma)][i(lambda)]
  op = list(outer.res=outer.res,gammas=gammas,lambdas=lambdas,n=length(eps.mimic), T=T,
            names=names,deltas=deltas,split=split,lambda.select.ix=lambda.select.ix,
            gamma.sel=gamma.sel,other.args=other.args)
  saveRDS(op,paste0(file,"op"))
  op
}

get.possible.selects = function(K,nDiverg){
  vecs = list()
  numbernonzero = 1
  combs=combn(K,nDiverg)
  for (c in 1:dim(combs)[2]){
    vec = rep(0,K) 
    for (i in 1:dim(combs)[1]){
      vec[combs[i,c]] = 1
    }
    vecs[[c]] = vec
  }
  vecs 
}

plot.lambda.paths = function(outer.res,gammaixleg=1,deltaixleg=3,max.n.div=1){
  #par(mfrow=c(2,1))
  
  ##### ADDED
  
  # plots lambda against coefficients
  other.args = outer.res$other.args

  max.n.div = other.args$C
  var.exists = !is.null(outer.res$or.shell.var)
  if (!var.exists){

    plot.emp.var = 0
  }else{
    outer.res.emp.var = outer.res$or.shell.var$outer.res
    plot.emp.var = other.args$plot.emp.var
    
  }
  ##### ADDED
  
  
  split=outer.res$split
  gamma.sel=other.args$gamma.select.ix
  delta.sel = other.args$delta.select.hardcode
  print(split)
# take c because sometimes it s a matrix if eg doing MC and averaging
  n.train = c(outer.res$outer.res[[1]][[1]][[1]]$n.train)
  n.test = c(outer.res$outer.res[[1]][[1]][[1]]$n.test)
  names = outer.res$names
  n = outer.res$n
  use.diff = outer.res$other.args$use.diff
  round.dig = 5
  if (split){
    n.train=n.train
    n.test=n.test
  }else{
    n.train=n
    n.test=n
  }
  deltas = outer.res$deltas
  print(n)
  lambda.select.ix = outer.res$lambda.select.ix
  lambdas=outer.res$lambdas
  gammas=outer.res$gammas
  outer.res= outer.res$outer.res
  nlam=length(lambdas)
  lambda.selects = list()
  jj=1
  #png("plts.png",width=500*8,height=500*2)
  value.ymin = NULL
  value.ymax = NULL
  
  roll.diff.p.gt.zs =  array(NA,c(length(deltas),length(gammas)))
  
  for (d in 1:length(deltas)){
    lambda.selects[[d]]=list()
    delta = deltas[d]
    dres = outer.res[[d]]
  
    ### ADDED
    
    if (var.exists){
      dres.emp.var = outer.res.emp.var[[d]]#
    }
    
    
    ### ADDED
    
  for (j in 1:length(gammas)){
    if (j==1){pick.lam=1}else{pick.lam=0}
    gamma = gammas[j]
    res = dres[[j]]
    betanlambdas.l=lapply(res,function(x){x[['betanlambda']]})
    betanlambdas=unlist.reshape(betanlambdas.l)
    betanlambdas.diff = betanlambdas
    
    meanpis = lapply(res,function(x){x[['rollout.betangammalambda.mean']]})
    meanpis = unlist.reshape(meanpis)

    diff.p.gt.z = lapply(res,function(x){x[['roll.diff.betangammalambda.bn.ab']]})
    diff.p.gt.z = unlist.reshape(diff.p.gt.z)
 
    mean.lp.rollout.betanlambda.sel = lapply(res,function(x){x[['mean.lp.rollout.betanlambda.sel']]})
    mean.lp.rollout.betanlambda.sel =  unlist.reshape(mean.lp.rollout.betanlambda.sel)
    
    nsel = lapply(res,function(x){x[['nsel']]})
    nsel =  unlist.reshape(nsel)

    mean.lp.rollout.betanlambda.nosel = lapply(res,function(x){x[['mean.lp.rollout.betanlambda.nosel']]})
    mean.lp.rollout.betanlambda.nosel = unlist.reshape(mean.lp.rollout.betanlambda.nosel)
    
    #mt=cbind(lambdas=lambdas,
    #      pibar=meanpis,
    #      probgtz=diff.p.gt.z,
    #      lp.sel=mean.lp.rollout.betanlambda.sel,
    #      lp.nosel=mean.lp.rollout.betanlambda.nosel)
    #browser()
    #bs.l = lapply(res)
    vars = unlist.reshape(lapply(res, function(x){diag(x$var.betanlambda)}))
    #eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T)
    #process.eps(eps,b0=NULL,gamma=1,u=1,delta=1)
    
    #must be vector for broadcasting later
    beh=c(res[[1]]$bn)
    K = length(beh)
 

    plot.adap.lasso.ci=0
    

    #zse = qnorm(0.975)*sqrt(vars/n.train)
    
    zse = sqrt(vars/n.train)

    #### ADDED
    
    if(var.exists){
      
      res.emp.var = dres.emp.var[[j]]
      
      vars.emp =  unlist.reshape(lapply(res.emp.var, function(x){x$betanlambda}))
      # we then construct ci for beta0 as betan +- sqrt(vars.emp)*z
      
      #zse.emp = qnorm(0.975)*sqrt(vars.emp)
      zse.emp = sqrt(vars.emp)

      mt=as.data.frame(cbind(emp=sqrt(vars.emp),theor=sqrt(vars/n.train)))
      colnames(mt) = c(paste0("empirical sd betanlambda",c(1:K)),
                       paste0("theoretical sd betanlambda",c(1:K)))
      rownames(mt)= lambdas
      print(mt)
      #browser()
      
    }
    
    ## ADDED
    
    plot.diff=0
    if (plot.diff){
      for (k in K:1){
        betanlambdas.diff[,k]=betanlambdas[,k]-beh[k]
        #zse[,k]=zse[,k]-beh[k]
      }
    }


    if (plot.diff){
      all=c(betanlambdas.diff)
    }else{
      all=c(beh,betanlambdas)
    }
    

    
    use.ci.for.ylim = other.args$use.ci.for.ylim
    if (use.ci.for.ylim){
    if (plot.adap.lasso.ci){
      
      for (i in 1:K){
      if (plot.diff){
        ylb=betanlambdas.diff[,i]-zse[,i]
        yub=betanlambdas.diff[,i]+zse[,i]
      }else{
        
       ylb = betanlambdas[,i]-zse[,i]
       yub = betanlambdas[,i]+zse[,i]
      }
        
        all = c(all,ylb,yub) 
      }
    }
    }
    logx = 0
    if(logx){
      logxarg = 'x'
    }else{
      logxarg = ''
    }

    
    #if(plot.diff){
    #  all = c()
    #  for (i in 1:K){
    #    all = c(all,betanlambdas[,i]-zse[,i],betanlambdas[,i]+zse[,i])
    #  }
    #}
    
    if (plot.diff){
      myylab=TeX("$\\beta_{n,\\gamma,\\lambda}-b_{n}$")
    }else{
      myylab=TeX("$\\beta_{n,\\gamma,\\lambda}$ or $b_{n}$")
    }
    
    # cofficient
    
    #bottom, left, top, and right
    par(mar=c(1,5.5,2,2))
    plot(lambdas,rep(0,nlam),ylim=c(min(all),max(all)+.2),lty=0,
         type='l',
         xaxt='n',
         cex.main=my.cex.main,
         cex.axis=my.cex.axis,
         cex.lab=my.cex.lab,
         #main=(bquote("n=" ~ .(n.train) ~ "; " ~ gamma == .(formatC(gamma, format = "e", digits = 2)) ~ '; ' ~ delta == .(delta)~ '; ' ~ Delta == .(use.diff))),
         main=(
           bquote(
             #"n=" ~ .(n.train) ~ "; " ~ 
                    gamma == .(round(gamma,2)
                                #formatC(gamma, 
                                       #format = "e",  #scientific notation
                                #       digits = 2)
                               ) ~ '; ' ~ delta == .(delta))),
         xlab=TeX("$\\lambda$"),
         ylab=myylab,log=logxarg)
    
    #text(min(lambdas),max(all)+0.2,"P(gtz)")
    #mtext('text is here', side=4, line=3.5, at=max(all))
    #text(lambdas,max(all)+.2,c(diff.p.gt.z))
    #text(min(lambdas),max(all)+.1,"nsel")
    #plot(lambdas,nsel,ylab=TeX("$d_{\\lambda_{n}}$"),xlab=TeX("$\\lambda$"))
    
   
    #mtext('text is here', side=4, line=3.5, at=max(all))
    #text(lambdas,max(all)+.1,c(nsel))
    
 
    
    colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    colorBlindBlack8 = palette.colors(palette = "Okabe-Ito")
    
    #### ADDED
    
    if (plot.emp.var){
      #for (i in 1:K){
      #  polygon(c(rev(lambdas),lambdas),c(rev(betanlambdas[,i]-zse.emp[,i]),
      #                                    betanlambdas[,i]+zse.emp[,i]),
      #          col=adjustcolor(colorBlindBlack8[i],alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
      #          lty=0)
      #}
      
      
      
      for (i in 1:K){
        if(plot.diff){
          beta.p = betanlambdas.diff[,i]
        }else{
          beta.p = betanlambdas[,i]
        }
        
        
        lines(lambdas,beta.p-zse.emp[,i],
                col=adjustcolor(colorBlindBlack8[i],alpha.f=1),#rgb(0,0,0,alpha=0.1),
                lty=3)
        lines(lambdas,beta.p+zse.emp[,i],
              col=adjustcolor(colorBlindBlack8[i],alpha.f=1),#rgb(0,0,0,alpha=0.1),
              lty=3)
      }
      
    }
    
    
    ## ADDED
    
    if (plot.adap.lasso.ci){
   
   # for (i in c(K:1)){
      if (var.exists){range=c(1:2)}else{range=c(1)}
      for (i in range){
      
      if (plot.diff){
        lb = betanlambdas.diff[,i]-zse[,i]
        ub =  betanlambdas.diff[,i]+zse[,i]
      }else{
      lb = betanlambdas[,i]-zse[,i]
      ub =  betanlambdas[,i]+zse[,i]
      }
      
      polygon(c(rev(lambdas),lambdas),c(rev(lb),
                                        ub),
              col=adjustcolor(colorBlindBlack8[i],alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
              lty=0)
    }
    }
    
   
    if (plot.diff){
    for(k in K:1){
      lines(lambdas,betanlambdas.diff[,k],lty=1,col=colorBlindBlack8[k],lwd=2)
    }
    }else{
    for (k in K:1){
      lines(lambdas,betanlambdas[,k],lty=2,col=colorBlindBlack8[k],lwd=2)
      abline(h=beh[k],lty=1,col=colorBlindBlack8[k],lwd=2)
    }
    }
    
    
    #vs = lapply(res,function(x){x$vn.betagammalambda})
    #vs = unlist(vs)
    
    #sigma.vs = lapply(res,function(x){x$sigam.vn.})
    #sigma.vs = unlist(sigma.vs)
    
    
    if(lambda.select.ix){
      #pick.lam =0 here avoids scale plot plotting lambda
      pick.lam=0
    }
    
    #lambda.star.ix=scale.plot(lambdas,sigma.vst,(vst[1]-vst),"test",plott=FALSE,pick.lam=pick.lam) #function(ls,x,y,label){
    #star.ix=lambda.star.ix$ix.star
    
    # if we set lambda select, use it, but overridden below
    if(lambda.select.ix){
      star.ix=lambda.select.ix
    }
    
    vs = lapply(res,function(x){x$vn.betagammalambda})
    ws = lapply(res,function(x){x$wn.betagammalambda})
    vs = unlist(vs)
    sigma.vs = lapply(res,function(x){x$sigam.vn.})
    sigma.vs = unlist(sigma.vs)
    #zse.v = qnorm(0.975)*sigma.vs/sqrt(n.train)
    zse.v = sigma.vs/sqrt(n.train)
    
    if (is.null(other.args$num.se)){
      
      vmin = other.args$vmin
    }else{
      # shoud just use behavioral and not bigest lambda.
      vmin = c(res[[length(res)]]$vn.betagammalambda + other.args$num.se*sqrt(res[[length(res)]]$sigam.vn./n.train))
    }
    
    
    
    
    select.lambda.pure.se = function(betanlambdas,beh,zse,use.diff,vs,vmin){
      # use.diff = Inf
      # # instead let's find lambda that satisfies our criterion
      # # of having a certain number diverging, but being the smallest lambda that does that
      if (use.diff){
        # diff big
        diff.from.behav = abs(t(betanlambdas)-beh)>use.diff
      }else{
        # diff big same as beh not in ci
        diff.from.behav = !(t(betanlambdas-zse)<=beh & t(betanlambdas+zse)>=beh)
      }
     
        
        gtvminix = which(vs>vmin)
        if (length(gtvminix)==0){
          print("nothing better than 2 se")
          gtvminix = c(1) #just choose first lambda (for debug mostly; in practice shouldn't happen)
        }
        
        star.ix = max(gtvminix)
        #star.ix = max(abs.min.gt.vmin)
        

      #}
      selected.mask=(diff.from.behav*1)[,star.ix]
      list(star.ix=star.ix,selected.mask=selected.mask)
    }
    
    select.lambda = function(betanlambdas,beh,zse,use.diff,vs,vmin){
      # instead let's find lambda that satisfies our criterion
      # of having a certain number diverging, but being the smallest lambda that does that
      if (use.diff){
        # diff big
        diff.from.behav = abs(t(betanlambdas)-beh)>use.diff
      }else{
        # diff big same as beh not in ci
        diff.from.behav = !(t(betanlambdas-zse)<=beh & t(betanlambdas+zse)>=beh)
      }
      select.lambda.adaptive = 1
      no.lambda.satisfies = 0
      if (select.lambda.adaptive){
        n.diff = apply(diff.from.behav,2,sum)
        # could possibly only take ones with positive here, so remove 000 selection
        # but in some ways 000 selection preferable "first do no harm
        
        dontLetEqualZero=0
        if (dontLetEqualZero){
          n.diff[n.diff==0] = 1e6
        }
        abs.diff.diff = abs(n.diff-max.n.div)
        
        #abs.min = which(abs.diff.diff==min(abs.diff.diff))
        
        
        gtvim  = vs>vmin
        # if it is not greater than minimum value, set it to inf
        abs.diff.diff[gtvim==FALSE] = Inf

        # take the min now of the ones, to get closest to target
        abs.min.gt.vmin = which((abs.diff.diff==min(abs.diff.diff)))
        
        # now take the max to get the most sparsity? 
        # or min to get most value. I think because
        # now we are over a value threshold, best to target sparsity
        # before, was no value threshold, so was taking largest value (min)
        # now we max. No now we still min
        star.ix = min(abs.min.gt.vmin)
        #star.ix = max(abs.min.gt.vmin)
        
        exact = 0
        if (exact){
          lam.sat.crit = (n.diff == max.n.div)
          if (sum(lam.sat.crit)==0){
            # if this is the case, no lambda satisfies our criterion. so set to 1
            # then i guess we just don't use these parameter choices ..
            # need to check for this...
            star.ix = length(lambdas)
            no.lambda.satisfies = 1
          }else{
            star.ix = min(which(lam.sat.crit==TRUE))
          }
        }
      }
      selected.mask=(diff.from.behav*1)[,star.ix]
      list(star.ix=star.ix,no.lambda.satisfies=no.lambda.satisfies,selected.mask=selected.mask)
    }
    

    sel.res = select.lambda.pure.se(betanlambdas,beh,zse,use.diff,vs,vmin=vmin) #set Vmin=-Inf to  just do coefficient
   
    star.ix = sel.res$star.ix
    no.lambda.satisfies = sel.res$no.lambda.satisfies
    selected.mask = sel.res$selected.mask
    
    # for now just use second delta, assuming always will have one larger and smallr
    if ((j==gamma.sel) & (d==delta.sel)){  
      abline(v=lambdas[star.ix],lty=2)
    }
 
    
    if ((j==gammaixleg) & (d==deltaixleg)){
   
    #if ((d==deltaixleg)){
    
      #legend("topright",paste0(rep(outer.res$names,each=2)," ",rep(c("b_n","\\beta_{n,\\gamma,\\lambda}"),each=2)),lty=rep(c(1,2),K),
    #col=rep(colorBlindBlack8[1:K],each=2),lwd=rep(2,2*K)) 
      
      #names= TeX(names)
      cols= colorBlindBlack8[1:K]
      lwd = rep(2,K)
      lty= rep(1,K)
      ncol=2
      text.width=1.2
      if(var.exists){
        names=c(names,"Emp.")
        cols=c(cols,1)
        lwd=c(rep(2,K),1)
        lty = c(lty,3)
        ncol=2
        text.width=NULL
      }
      #$\\beta_{n,\\gamma,\\lambda}$ or $b_{n}$
      
      names=c(c("$\\beta_{n,\\gamma,\\lambda}$","$b_{n}$"),names)
      lwd=c(c(2,2),lwd)
      lty=c(c(2,1),lty)
      cols = c(c(1,1),cols)
     print("printing Legend") 
    legend("topright",TeX(names),lty=lty, ncol=ncol,#inset=c(-0.2,0),
           col=cols,lwd=lwd,bg="white",cex=cex.legend)#,text.width=1.2)
    }
    


    
    ##### ADDED
    
    if(var.exists){
      
      res.emp.var = dres.emp.var[[j]]

      M = length(res.emp.var)


      vars.v.emp =  unlist.reshape(lapply(res.emp.var, function(x){x$vn.betagammalambda}))
      # we then construct ci for beta0 as betan +- sqrt(sigma^2/n)*z
      #zse.v.emp = qnorm(0.975)*sqrt(vars.v.emp)
      zse.v.emp = sqrt(vars.v.emp)

      mt=as.data.frame(cbind(emp=sqrt(vars.v.emp),theor=sigma.vs/sqrt(n.train)))
      colnames(mt) = c("empirical sd Vn","theoretical sd Vn")
      rownames(mt)= lambdas
      print(mt)
      #browser()

      
          }

    
    ##### ADDED
    
    
    
    ## test
    
    vst = lapply(res,function(x){x$vn.betagammalambda.test})
    wst =  lapply(res,function(x){x$wn.betagammalambda.test})
    vst = unlist(vst)
    sigma.vst = lapply(res,function(x){x$sigam.vn.test})
    sigma.vst = unlist(sigma.vst)
    #zse.vt = qnorm(0.975)*sigma.vst/sqrt(n.test)
    zse.vt = sigma.vst/sqrt(n.test)
    
    #plt.train=1
    #if (plt.train){
    plot.w = 0
    if (plot.w){
    ylab = TeX('$V_{n}$ or $W_{n}$')
    }else{
      ylab = TeX('$V_{n}$')
    }
    use.first.val.for.ylim = 0
    use.value.se.for.ylim = other.args$use.value.ci.for.ylim
    #use.value.se.for.ylim = 0#1
    if (use.first.val.for.ylim){
    if (j==1 & d ==1){
      # take the value bounds from the first, and use them for every plot
      
      if (use.value.se.for.ylim){
      value.ymax = max(c(vs+zse.v,vst+zse.vt,vs))
      value.ymin = min(c(vs-zse.v,vst-zse.vt,vs))
      }else{
        value.ymax = max(c(vs,vst))
        value.ymin = min(c(vs,vst))
        
      }
      
    }
    }else{
      if (use.value.se.for.ylim){
        value.ymax = max(c(vs+zse.v,vst+zse.vt,vs))
        value.ymin = min(c(vs-zse.v,vst-zse.vt,vs))
      }else{
        value.ymax = max(c(vs,vst))
        value.ymin = min(c(vs,vst))
        
      }
      
    }
    
    
    # diff pi mean
    
    #bottom, left, top, and right
    par(mar=c(0.5,5.1,1,2))
    
    plot(lambdas,diff.p.gt.z, 
         #ylab=TeX("$| \\bar{\\pi}_{1}(\\beta_{n,\\gamma,\\lambda}) - \\bar{\\pi}_{1}(b_{n})|$"),
         ylab=TeX("$|\\bar{\\pi}_{sugg}-\\bar{\\pi}_{beh}|$"),
         #ylab=TeX("$| P_{sugg}(Tx) - P_{beh}(Tx)$"),
         xlab=TeX("$\\lambda$"),
         main="",
         xaxt='n',
         #main="Mean probability treatment",
         log=logxarg,type='l',
         cex.main=my.cex.main,
         cex.axis=my.cex.axis,
         cex.lab=my.cex.lab)
    
    if ((j==gamma.sel) & (d==delta.sel)){  
      abline(v=lambdas[star.ix],lty=2)
    }
    
    #### value
  
    
    #bottom, left, top, and right
    par(mar=c(4,5.1,2,2))
    
    plot(lambdas,vs,ylab=ylab,xlab=TeX("$\\lambda$"),type='l',lwd=3,
           ylim=c(value.ymin,value.ymax),
         log=logxarg,
         main="",
         cex.main=my.cex.main,
         cex.axis=my.cex.axis,
         cex.lab=my.cex.lab,
  
         #main="Value",
         #main=bquote("n=" ~ .(n.test))
         #main=(bquote("n=" ~ .(n.train) ~ "; " ~ gamma == .(formatC(gamma, 
#format = "e", digits = 2)) ~ '; ' ~ delta == .(delta)~ '; ' ~ Delta == .(use.diff))),
         )

    if (plot.w){
    lines(lambdas,ws,lty=2,col=3)
    }
    #lines(lambdas,vs-zse.v,lty=2)

    polygon(c(rev(lambdas),lambdas),c(rev(vs-zse.v),
                                        vs+zse.v),
              col=adjustcolor(1,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
              lty=0)
    
    ### ADDED
    if (var.exists){
    lines(lambdas,vs-zse.v.emp,lty=3)
    lines(lambdas,vs+zse.v.emp,lty=3)
    }
    ### ADDED
    
    #lower.bound = vs-zse.v
    #star.ix = which(lower.bound==max(lower.bound))
    #points(lambdas[star.ix],lower.bound[star.ix],pch=8)
    
    #star.ix.v = which(vs==max(vs))
    #points(lambdas[star.ix],vs[star.ix],pch=8)
 
    
    #scale.plot(lambdas,sigma.vs,(vs[1]-vs),"train",plott=FALSE) #function(ls,x,y,label){

    lines(lambdas,vst,lty=3,col=2,lwd=3)#,ylab=TeX('$V_{n}(test)$'),xlab=TeX("$\\lambda$"),type='l',
         #ylim=c(min(c(vs-zse.v,vs)),max(c(vs+zse.v,vs))))
    if (plot.w){
    lines(lambdas,wst,lty=4,col=4)
    }

    #lines(lambdas,vs-zse.v,lty=2)
    polygon(c(rev(lambdas),lambdas),c(rev(vst-zse.vt),
                                      vst+zse.vt),
            col=adjustcolor(2,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
            lty=0)
    lower.bound = vst-zse.vt
    
    if ((j==gammaixleg) & (d==deltaixleg)){
    #if ((d==deltaixleg)){
      if (plot.w){
        legend("topright",TeX(c("Train $V_n$","Test $V_n$","Train $W_n$", "Test $W_n$")),
               lty=c(1,3,2,4),col=c(1,2,3,4),lwd=c(3,3,1,1))
      }else{
        
        names = c("Train","Test")
        lty = c(1,3)
        col = c(1,2)
        lwd=c(2,2)
        ncol=2
        if (var.exists){
          names=c(names,"Emp.")
          lty = c(lty,3)
          col=c(col,1)
          ncol=3
          lwd=c(lwd,1)
          
          
        }
        
        legend("topright",names,lty=lty,col=col,ncol=ncol,
               bg="white",lwd=lwd,cex=cex.legend,y.intersp=0.5,x.intersp=0.5,text.width=0.1)
        only.micu=1
        if(only.micu){
          legend("topright",names,lty=lty,col=col,ncol=ncol,
                 bg="white",lwd=lwd,cex=cex.legend)
        }
        
      }
    }
    
    #star.ix = which(lower.bound==max(lower.bound))
    #points(lambdas[star.ix],lower.bound[star.ix],pch=8)
    #text(lambdas[star.ix],lower.bound[star.ix],"V-2sigma/sqrt(n)")
    
    #star.ix.v = which(vst==max(vst))
    #points(lambdas[star.ix],vs[star.ix],pch=8)
   # text(lambdas[star.ix],vs[star.ix],"V",offset=.01)
    #plot(lambdas,sigma.vs,ylab=TeX('$\\sigma(V_{n}(test (NOT TEST)))$'),xlab=TeX("$\\lambda$"),type='l')  
    
   
    # for now just use second delta, assuming always will have one larger and smallr
    if ((j==gamma.sel) & (d==delta.sel)){  
      abline(v=lambdas[star.ix],lty=2)
    }
    
    # take the first lambdas roll diff
    # sometimes though coefficients can increase with increasing lamdba
    # take lambda=1
    
    #
    smallest.lambda.ix = 1
    roll.diff.p.gt.zs = res[[smallest.lambda.ix]]$roll.diff.p.gt.z
    
    print(c("gamma",gammas[j],"delta",deltas[d]))

    print(c("ix selected:",star.ix))
    print(c("selected mask:",selected.mask))

    lambda.selects[[d]][[j]]=list(lambda.star.ix=star.ix,
                              delta=delta,gamma=gamma,
                              selected.mask=selected.mask, 
                              no.lambda.satisfies=no.lambda.satisfies,
                              roll.diff.p.gt.zs=roll.diff.p.gt.zs)
    jj=jj+1
  }
  }
  lambda.selects
}

select.and.est = function(eps,b0,gammas,lambdas,deltas,names,resfile,plotfile,scale.s,
                          lambda.select.ix=FALSE,gamma.select.ix=2,rdigits=3,other.args){
  K=dim(eps[[1]]$Ss[[1]])[1]
 
  sel = rep(1,K)
  
  n = length(eps)
  n.train = as.integer(length(eps)/2)
  n.test = n-n.train
  #n.test = outer.res$outer.res[[1]][[1]][[1]]$n.test
  #n.train =  outer.res$outer.res[[1]][[1]][[1]]$n.train
  #n=outer.res$outer.res[[1]][[1]][[1]]$n.test
  train.eps = eps[1:n.train]
  test.eps=eps[(n.train+1):n]

  # should in some ways move split out of this function
  outer.res = lambda.path(train.eps,b0=b0,gammas=gammas,lambdas=lambdas,
                          names=names,
                          #names=c(TeX('$S_1$'),TeX('$S_2$')),
                          deltas,file=paste0(resfile,"mimic.outer.res"),sel=sel,split=1,scale.s=scale.s,
                          lambda.select.ix=lambda.select.ix,gamma.sel=gamma.select.ix,
                          other.args=other.args)
  
  print(paste0("rollout of selected policy,
               outer res is outer.res[d(delta)][j(gamma)][i(lambda)]"))
  print("note that lambda does change here")
  print(outer.res$outer.res)
  
  
  #mc.res[[i]]=outer.res
  #}
  ngam = length(gammas)
  nlam = length(lambdas)
  ndel = length(deltas)
                
  #outer.ress=lapply(mc.res,function(x){x$outer.res})
  #lapply(outer.ress,function(x){x[[2]]})
  #png("grid.mc.png",height=400*7*ngam,width=500*ndel,res=200)
  #saveRDS(outer.res,resfile)
  #outer.res = readRDS(resfile)
  
  #}

  nplots = 3
  
  res = ngam*res.multiplier

  png(plotfile,height=heightchunk*ngam*nplots,width=widthchunk*ndel,res=res)
  #png(plotfile,height=500*ngam*nplots,width=900*ndel,res=res)
  #par(mfcol=c(ngam*nplots,ndel)) #vertical in gamma
  
  

  ml=layout(matrix(c(1:(ndel*nplots*ngam)),ncol=ndel), 
            widths=rep(myw,nplots*ngam*ndel), 
            heights=rep(myheights,ngam), TRUE)
  

  
  #par(mfcol=c(ngam*5,ndel)) #horiz in gamma
  # : bottom, left, top, and right.
  
  ixs=plot.lambda.paths(outer.res,max.n.div=other.args$C)
  dev.off()
  
  mar=get.stats(outer.res)
  saveRDS(mar,paste0(resfile,'stats'))
  
  ratios  = get.ratios(outer.res)
  saveRDS(ratios,paste0(resfile,'ratios'))
  
  #if (lambda.select.ix){
  #  select.mask = lambda.select.ix
  #}else{
  
  # hard coding delta=2 selection
  delta.select.hardcode = other.args$delta.select.hardcode
  gamma.select.ix = other.args$gamma.select.ix
  #delta.select.hardcode= 1#2
  select.mask =  ixs[[delta.select.hardcode]][[gamma.select.ix]]$selected.mask
  #}
  
  # or is [[delta]][[gamma]][[lambda]]
  #n.test = outer.res$outer.res[[1]][[1]][[1]]$n.test
  #n.train =  outer.res$outer.res[[1]][[1]][[1]]$n.train
  #test.eps=eps[(n.train+1):length(eps)]
  
  outer.res.sel = lambda.path(test.eps,b0=b0,gammas=gammas,lambdas=c(0),
                              names=names,
                              #names=c(TeX('$S_1$'),TeX('$S_2$')),
                              deltas,file=paste0(resfile,"mimic.outer.res.sel"),
                              sel=select.mask,split=0,scale.s=scale.s,
                              lambda.select.ix=lambda.select.ix,
                              gamma.sel=gamma.select.ix,
                              other.args=other.args)

  n.sel=outer.res.sel$n
  print(c("correct n (check)?",n.sel==n.test))
  # we only take half the data
 
  # selection ix corresponds to which gamma we select for
  #selection.ix = gamma.select.ix
  
  # it's outer.res[[delta]][[gamma]][[lambda]]. lambda is 1, which is 0, since we are doing this
  # after selection.
  
  #for (d in 1:length(deltas)){
  # does not depend on delta

  d=1
  for (selection.ix in 1:length(gammas)){
  
  cbs=sqrt(diag(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$var.betanlambda))/sqrt(n.sel)*qnorm(0.975)
  cbs.b = sqrt(diag(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$var.b))/sqrt(n.sel)*qnorm(0.975)
  ests= c(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$betan)
  ests.b = c(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$bn)
  

  coef.table = cbind(beta=ests,lcb=(ests-cbs),ucb=(ests+cbs),sel=select.mask,
                     bn=outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$bn)
  rownames(coef.table)= outer.res$names#paste0("$S_",1:K,"$")
  colnames(coef.table)=c("est", "LCI","UCI","Active","$b_n$")
  rdigits = rdigits
  coef.table=formatC(coef.table,format='f',digits=rdigits)

  coef.table=as.data.frame(coef.table)
  beta.and.ci=paste0(coef.table$est," (",coef.table$LCI,", ",coef.table$UCI,")")
  b.and.ci = paste0(formatC(ests.b,format="f",digits=rdigits),
                    " (",formatC(ests.b-cbs.b,format='f',digits=rdigits),
                    ", ",formatC(ests.b+cbs.b,format="f",digits=rdigits),")")
  newcoef.table=cbind(beta.and.ci,b.and.ci,select.mask)
  # I think should be ok
  colnames(newcoef.table)=c("Suggested ($\\beta_{n,\\gamma}$)","Behavioral ($b_n$)","Active")
  rownames(newcoef.table)= outer.res$names
  print(newcoef.table)
  newcoef.table = as.data.frame(newcoef.table)

  #print(xtable(newcoef.table[newcoef.table[,'Active']==1,1:2], caption = paste0("$n=$",n.sel,", $\\gamma=$",
  newcoef.table.set = newcoef.table
  newcoef.table.set[newcoef.table.set[,'Active']==0,1] = "set to $b_n$"
  # I hope new caption is ok
  print(xtable(newcoef.table.set[,1:2], caption = paste0("Estimated coefficients (95$\\%$ CI) for held out dataset post-selection inference; $n=$",n.sel,", $\\gamma=$",
                                               #round(gammas[selection.ix],4)
                                               #formatC(gammas[selection.ix], format = "e", digits = 2), ", $\\delta=$",deltas[d])),
                                               formatC(gammas[selection.ix], format = "e", digits = 2))),
        type="latex",sanitize.text.function = function(x){x},include.rownames=TRUE,digits=4)
  }
 # }
  
  print(paste0("rollout of selected policy,
               outer res is outer.res[d(delta)][j(gamma)][i(lambda)]"))
  print("note that lambda doesnt change here")
  print(outer.res.sel$outer.res)
  # note that n comes from now sel, which should be just half data

  list(ixs=ixs,outer.res=outer.res)
}

check.pos = function(eps.mimic,cov.of.int){
  T = length(eps.mimic[[1]]$As)
  sc=center.scale.s(eps.mimic)
  a=unlist(lapply(sc$e,function(x){x$As}))
  r=unlist(lapply(eps.mimic,function(x){x$Rs}))
  smm=sc$sm.cs
  #apply(smm,2,var)
  
  colnames(smm)=cov.of.int
  beh = glm(a~smm-1,family='binomial')
  pi = beh$fitted.values
  n = length(eps.mimic)
  na.states=apply(is.na(smm),1,sum)
  
  #knn(train_set=as.data.frame(cbind(a=a,smm),test_set=as.data.frame(smm),k=2,categorical_target = "a"))
  df.d=as.data.frame(cbind(r,a,smm))
  par.te=summary(lm(r~.,data=df.d))
  mu.1 = 1/n*sum(a*r/(pi))
  mu.0 = 1/n*sum((1-a)*r/(1-pi))
  treat.eff = mu.1 - mu.0
  df = as.data.frame(cbind(obs=a,pred=pi))
  #calibration_plot(data=df,obs="a",pred="pred")
  prop.r.less.than.0=mean(r<0)
  prop.a.1 = mean(a)
  #print(summary(beh))
  #print(xtable(summary(beh)))
  #pi
  #hist(pi)
  #hist(pi,main=careunit)
  png("Overlap.png",width=1000,height=dim(msr)[2]*1000,res=150)
  min.beh=min(c(1-pi,pi))
  max.beh=max(c(1-pi,pi))
  mean.beh=mean(c(1-pi,pi))
  hist(pi[as.logical(1-a)],col=rgb(.9,.9,.9),xlab=TeX("$\\pi_{b_n}(A_0=1|s)$"),main="Overlap")
  hist(pi[as.logical(a)],add=TRUE,col=rgb(.5,.5,.5))
  legend("topright",lty=c(1,1),lwd=c(4,4),c("Observed A=0","Observed A=1"),col=c(rgb(.9,.9,.9),rgb(.5,.5,.5)))
  dev.off()
  tb=c("n"=as.integer(n),
       "T"=as.integer(T),
       "Prop. R<0"=round(prop.r.less.than.0,3),
       "Prop A=1"=round(prop.a.1,3),
       "min $\\pi_b$"=round(min.beh,3),
       "mean $\\pi_b$"=round(mean.beh,3),
       "max $\\pi_b$"=round(max.beh,3),
       "IPTW treat. eff."=round(treat.eff,3),
       "Par. treat. eff (p-val)"=paste0(round(par.te$coefficients['a','Estimate'],3),
                                        " (",round(par.te$coefficients['a','Pr(>|t|)'],3),")")
  )
  
  dftb=as.data.frame(tb)
  print(xtable(dftb),
        type="latex",sanitize.text.function = function(x){x},
        include.rownames=TRUE,digits=2)
  
  list(n=length(eps.mimic),T=T,
       min.beh=min.beh,
       max.beh=max.beh,
       mean.beh=mean.beh,
       summ=summary(beh),IPTW.treat.eff=treat.eff,
       par.treat.eff=par.te$coefficients['a',],
       prop.r.less.than.0=prop.r.less.than.0,
       prop.a.1=prop.a.1,
       dftb=dftb,beh=beh)
}


check.pos.cal = function(eps.mimic,cov.of.int,pen.b,train.test=1,plot.=TRUE,usesampler=FALSE,other.args){
  # checks positivity, etc for behavioral policy in real data analysis
  
  
  n = length(eps.mimic)
  
  T = length(eps.mimic[[1]]$As)
  
  
  # split data into test and training
  # if(split){
  # need to deal with this
  
  hf = as.integer(n/2)
  
  browser()
  #stage.randoms = rep(NA,)
  
  
  
  
  if (train.test){
    
    if (usesampler){  
      #train.ixs = sample(1:n,hf)  
      train.ix = as.logical(rbinom(n,1,prob=0.5))
      eps.test = eps.mimic[!train.ix]#[(hf+1):n]
      eps.train = eps.mimic[train.ix]#[1:hf]
    }else{
      eps.test = eps.mimic[(hf+1):n]
      eps.train = eps.mimic[1:hf]
      
    }
    
  }else{
    eps.test = eps.mimic
    eps.train = eps.mimic
  }
  
  
  n.test = length(eps.test)
  n.train = length(eps.train)
  #  }else{
  #    eps.test=eps
  #  }
  
  # if(scale.s){
  so.train = center.scale.s(eps.train)
  eps.train = so.train$e
  
  so.test = center.scale.s(eps.test,means=so.train$means,sds=so.train$sds)
  eps.test = so.test$e
  
  # }else{
  #   # don t center and scale
  #   so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K)) #not centering and scaling
  #   # but we still need the sds for the scaling factors.
  #   soc = center.scale.s(eps)
  #   sigmas = soc$sds
  #   # generally, we would do the above with scale factors. instead
  #   # override scale factors by setting sigmas to 1. we are not scaling for simulations. 
  #   # we just generate data with variance 1
  #   #we just scale for real data
  #   sigmas = rep(1,K)
  # }
  
  
  # first scale episodes
  #sc=center.scale.s(eps.mimic)
  
  get.matrices  = function(eps,sc,cov.of.int){  
    a=unlist(lapply(eps,function(x){x$As}))
    r=unlist(lapply(eps,function(x){x$Rs}))
    smm=sc$sm.cs[,1:length(cov.of.int),drop=FALSE]
    colnames(smm)=cov.of.int
    list(a=a,r=r,smm=smm)
  }
  
  d.train = get.matrices(eps.train,so.train,cov.of.int)
  a.train = d.train$a
  r.train = d.train$r
  smm.train = d.train$smm 
  
  d.test = get.matrices(eps.test,so.test,cov.of.int)
  a.test = d.test$a
  r.test = d.test$r
  smm.test = d.test$smm 
  
  
  #hf = 
  
  if(pen.b){
    
    K = dim(smm.train)[2]
    #lambdas = seq(0,.01,length.out=100)
    
    beh.train=cv.glmnet(smm.train,a.train,intercept = FALSE,alpha=0,
                        family='binomial',
                        type.measure=other.args$type.measure,
                        lambda=other.args$b.lambdas)
    
    
    
    beh.coef=coef(beh.train,s = "lambda.min")[2:(K+1),]
    
    #beh.coef=coef(beh.train,s = "lambda.1se")[2:(K+1),]
    
    
    #print(beh.coef)
    pi.train = expit(smm.train%*%beh.coef)
    pi.test = expit(smm.test%*%beh.coef)
    
    
  }else{
    beh.train = glm(a.train~smm.train-1,family='binomial')
    
    pi.train = expit(smm.train%*%beh.train$coefficients)#beh.train$fitted.values
    pi.test = expit(smm.test%*%beh.train$coefficients)
  }
  
  
  
  #dev.off()
  
  
  
  get.treatment= function(r,a,smm,pi){
    
    n = length(eps.mimic)
    #na.states=apply(is.na(smm),1,sum)
    #knn(train_set=as.data.frame(cbind(a=a,smm),test_set=as.data.frame(smm),k=2,categorical_target = "a"))
    df.d=as.data.frame(cbind(r,a,smm))
    par.te=summary(lm(r~.,data=df.d))
    mu.1 = 1/n*sum(a*r/(pi)) 
    mu.0 = 1/n*sum((1-a)*r/(1-pi))
    treat.eff = mu.1 - mu.0
    df = as.data.frame(cbind(obs=a,pred=pi))
    #calibration_plot(data=df,obs="a",pred="pred")
    prop.r.less.than.0=mean(r<0)
    prop.a.1 = mean(a)
    
    # overlap plot
    #png("Overlap.png",width=1000,height=dim(msr)[2]*1000,res=150)
    #min.beh=min(c(1-pi,pi))
    #max.beh=max(c(1-pi,pi))
    #mean.beh=mean(c(1-pi,pi))
    min.beh=min(c(1-pi,pi))
    max.beh=max(c(1-pi,pi))
    mean.beh=mean(pi)
    
    
    list(par.te=par.te,treat.eff=treat.eff,prop.r.less.than.0=prop.r.less.than.0,
         prop.a.1=prop.a.1,min.beh=min.beh,max.beh=max.beh,mean.beh=mean.beh)
  }
  
  d.train=get.treatment(r.train,a.train,smm.train,pi.train)
  par.te.train = d.train$par.te
  treat.eff.train = d.train$treat.eff
  
  d.test=get.treatment(r.test,a.test,smm.test,pi.test)
  par.te.test = d.test$par.te
  treat.eff.test = d.test$treat.eff
  
  
  
  
  prop.r.less.than.0.train=d.train$prop.r.less.than.0
  prop.a.1.train = d.train$prop.a.1
  min.beh.train = d.train$min.beh
  max.beh.train = d.train$max.beh
  mean.beh.train = d.train$mean.beh
  
  
  
  #par(mfrow=c(2,1))
  #png(paste0("pen.b=",pen.b,"CalibrationTrain.png"),width=500,height=500,res=100)
  my.data.train = data.frame(a.train,pi.train)
  p.train=calibration_plot(data=my.data.train,obs="a.train",pred="pi.train",nTiles=5,title="Train Calibration",data_summary = TRUE)
  if (plot.){
    print(p.train)
    #dev.off()
  }
  #png(paste0("pen.b=",pen.b,"CalibrationTest.png"),width=500,height=500,res=100)
  my.data.test = data.frame(a.test,pi.test)
  p.test=calibration_plot(data=my.data.test,obs="a.test",pred="pi.test",nTiles=5,title="Test Calibration",data_summary = TRUE)
  if (plot.){
    print(p.test)
  }
  
  if (plot.){
    hist(pi.train[as.logical(1-a.train)],col=rgb(.9,.9,.9),xlab=TeX("$\\pi_{b_n}(A_0=1|s)$"),main="Overlap (Training set)")
    hist(pi.train[as.logical(a.train)],add=TRUE,col=rgb(.5,.5,.5))
    legend("topright",lty=c(1,1),lwd=c(4,4),c("Observed A=0","Observed A=1"),col=c(rgb(.9,.9,.9),rgb(.5,.5,.5)))
    #dev.off()
  }
  
  
  
  tb=c("n"=as.integer(n),
       #"T"=as.integer(T),
       "Prop. R<0"=round(prop.r.less.than.0.train,3),
       "Prop A=1"=round(prop.a.1.train,3),
       "min $\\pi_b$"=round(min.beh.train,3),
       "mean $\\pi_b$"=round(mean.beh.train,3),
       "max $\\pi_b$"=round(max.beh.train,3),
       "IPTW treat. eff. train"=round(treat.eff.train,3),
       "IPTW.treat.eff.test"=round(treat.eff.test,3),
       "Par treat eff train (pval)"=paste0(round(par.te.train$coefficients['a','Estimate'],3)," (",round(par.te.train$coefficients['a','Pr(>|t|)'],3),")"),
       "Par. treat. eff test (p-val)"=paste0(round(par.te.test$coefficients['a','Estimate'],3)," (",round(par.te.test$coefficients['a','Pr(>|t|)'],3),")")
  )
  
  dftb=as.data.frame(tb)
  #print(xtable(dftb),
  #      type="latex",sanitize.text.function = function(x){x},
  #      include.rownames=TRUE,digits=2)
  
  list(n=length(eps.mimic),
       #T=T,
       min.beh.train=min.beh.train,
       max.beh.train=max.beh.train,
       pi.train=pi.train,
       pi.test=pi.test,
       a.train=a.train,
       a.test=a.test,
       mean.beh.train=mean.beh.train,
       summ=summary(beh.train),
       IPTW.treat.eff.train=treat.eff.train,
       par.treat.eff.test=par.te.train$coefficients['a',],
       IPTW.treat.eff.test=treat.eff.test,
       par.treat.eff.test=par.te.test$coefficients['a',],
       prop.r.less.than.0.train=prop.r.less.than.0.train,
       prop.a.1.train=prop.a.1.train,
       dftb=dftb,beh.train=beh.train,
       p.train=p.train$data_summary,
       p.test=p.test$data_summary)
}


plot.objs = function(eps,gammas,u,b0){
  T = length(eps[[1]]$As)
  ngam=length(gammas)
  par(mfrow=c(ngam,1))
  for (i in 1:ngam){
  #plot.d.mdp(eps)
  #title("scaled states")
  K=1
  #par(mfrow=c(1,1))
  as=unlist.reshape(lapply(eps,function(x){x$As}))
  mns=apply(as,2,mean)
  beta.grid = seq(-10,10,length.out=50)
  vn.pds = rep(NA,length(beta.grid))
  mn.pds = rep(NA,length(beta.grid))
  vns = rep(NA,length(beta.grid))
  mns = rep(NA,length(beta.grid))
  mns2 = rep(NA,length(beta.grid))
  v0s = rep(NA,length(beta.grid))
  v0s.gamma0 =  rep(NA,length(beta.grid))
  
  #v0s.gamma1e6 =  rep(NA,length(beta.grid))
  for (j in 1:length(beta.grid)){
    # would be interesting to plot wn.
    mns[j]=mn(beta=beta.grid[j],b=b0,eps=eps,
              gamma=gammas[i],u=u,sigmas=rep(1,K),sel=rep(1,K))
    #mns2[j]=mn(beta=beta.grid[j],b=b0,eps=eps,
    #           gamma=gamma,u=u,sigmas=so$sds,sel=rep(1,K))
    #mn.pds[j] = mn.pd(beta=beta.grid[j],b=b0,eps=eps.sc,gamma=gamma,u=u,
    #                  sigmas=rep(1,K),sel=rep(1,K))
    ###
    #vn.pds[j]=vn.pd(beta=beta.grid[j],b=b0,eps=eps,u=u)
    ###
    #vns[j] = vn(beta=beta.grid[j],b=b0,eps=eps,u=u,sel=rep(1,K))
    ###
    v0s.gamma0[j]=est.m0.lambda(beta=beta.grid[j],eps=eps,
                                gamma=gammas[i],u=u,b=b0,sigmas=rep(1,K),
                                sel=rep(1,K))
    
    # it's just an upside down cup centered at 0, much deeper than others
    #v0s.gamma1e6[j]=est.m0.lambda(beta.grid[j],eps,gamma=1e6,u) 
  }

  at=c(mns,v0s.gamma0)#,v0s.gamma1e6)
  #print("OK")
  #plot(beta.grid,mns)
  #lines(beta.grid*so$sds,mns2)
  
  tcol = rgb(0,0,0,alpha=0.5)
  ecol = rgb(0,1,1,alpha=0.5)
  ecol2 = rgb(0,0,1,alpha=0.5)
  plot(beta.grid,v0s.gamma0,lty=1,type='l',col=tcol,
       ylim=c(min(at),max(at)),xlab=expression(beta),ylab="M",
       main=paste0("gamma=",gammas[i],"; n=",length(eps),"T=",T))
  lines(beta.grid,mns,lty=2,lwd=4,col=ecol)
  #lines(beta.grid,mn.pds,lty=5,lwd=1,col=ecol2)
  
  #lines(beta.grid,v0s.gamma1e6,lty=4,lwd=0.5,col=rgb(0,1,1,alpha=1))
  if (i == 1){
  legend("topright",
         c(expression("M_0","M_n","M_n(PD)")),
         lty=c(1,2,5),
         lwd=c(1,4,1),
         col=c(tcol,
               ecol,
               ecol2))
  }
  #all = c(v0s.gamma0,vns)
  #plot(beta.grid,v0s.gamma0,lty=3,lwd=2,col=tcol,type='l',ylim=c(min(all),max(all)))
  #lines(beta.grid,vn.pds,lty=2,lwd=2,col=tcol)
  #lines(beta.grid,vns,lty=4,lwd=2,col=tcol)
  #legend("bottomleft",
  #       c(expression("V0"),"Vn(PD)","Vn"),#,expression(gamma~"=1e6")),
  #       lty=c(3,2,4),#,4),
  #       lwd=c(2,2,2),#,.5),
  #       col=c(tcol,
  #             tcol,
  #             tcol
               #rgb(0,1,1,alpha=1)
  #       ))
  }
  
}

get.eps.ms = function(s,a,r){
  
  # s is ldf
  n=dim(s[[1]])[1]
  a = as.matrix(a)
  r = as.matrix(r)
  T = length(a[1,])
  print(c("T",T))
  eps=list()
  
  for (i in 1:n){
    Ss=list()
    As=list()
    Rs=list()  
    for (t in 1:T){
      Ss[[t]]=matrix(as.vector(s[[t]][i,])) # why
      As[[t]]=a[i,t]
      Rs[[t]]=r[i,t] 
    }
    eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
}

get.eps = function(s,a,r){
  n=dim(s)[1]
  T = length(a[1])
  eps=list()
  
  for (i in 1:n){
    Ss=list()
    As=list()
    Rs=list()  
    for (t in 1:T){
      Ss[[t]]=matrix(s[i,]) # why
      As[[t]]=a[i]
      Rs[[t]]=r[i] 
    }
    eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
}


# plot.lambda.paths.MC=function(outer.res,or.shell.var=NULL){
#   
#   #par(mfrow=c(2,1))
#   names = outer.res$names
#   n = outer.res$n
#   deltas = outer.res$deltas
#   print(n)
#   lambdas=outer.res$lambdas
#   gammas=outer.res$gammas
#   outer.res= outer.res$outer.res
#   if(!is.null(or.shell.var)){
#     outer.res.var = or.shell.var$outer.res
#   }
#   nlam=length(lambdas)
#   #png("plts.png",width=500*8,height=500*2)
#   for (d in 1:length(deltas)){
#     delta = deltas[d]
#     dres = outer.res[[d]]
#     if(!is.null(or.shell.var)){
#       dres.var = outer.res.var[[d]]
#     }
#     for (j in 1:length(gammas)){
#       gamma = gammas[j]
#       res = dres[[j]]
#       if(!is.null(or.shell.var)){
#         res.var = dres.var[[j]]
#       }
#       bnls.l=lapply(res,function(x){x[['betanlambda']]})
#       bnls=unlist.reshape(bnls.l)
#       #bs.l = lapply(res)
#       vars = unlist.reshape(lapply(res, function(x){diag(x$var.betanlambda)}))
#       if(!is.null(or.shell.var)){
#         vars.n = unlist.reshape(lapply(res.var, function(x){(x$betanlambda)}))
#       }
#       
#       #eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T)
#       #process.eps(eps,b0=NULL,gamma=1,u=1,delta=1)
#       beh=res[[1]]$bn
#       all=c(beh,bnls)
#       plot(lambdas,rep(0,nlam),ylim=c(min(all),max(all)),lty=0,
#            type='l',main=(bquote(gamma == .(gamma) ~ '; ' ~ delta == .(delta))),
#            xlab=TeX("$\\lambda$"),
#            ylab=TeX("$\\beta_{n,\\gamma,\\lambda}$ or $b_{n}$"))
#       K = length(res[[1]]$bn)
#       
#       zse = qnorm(0.975)*sqrt(vars/n)
#       if(!is.null(or.shell.var)){
#         zse.n = qnorm(0.975)*sqrt(vars.n/M)
#       }
#       for (i in 1:K){
#         polygon(c(rev(lambdas),lambdas),c(rev(bnls[,i]-zse[,i]),
#                                           bnls[,i]+zse[,i]),
#                 col=adjustcolor(i,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
#                 lty=0)
#         if(!is.null(or.shell.var)){
#           # somehow wrong dim
#           
#           lines(lambdas,bnls[,i]-zse.n[,i],lty=2,col=3)
#           lines(lambdas,bnls[,i]+zse.n[,i],lty=2,col=3)
#         }
#       }
#       
#       for (k in K:1){
#         lines(lambdas,bnls[,k],lty=k,col=k,lwd=1.5)
#         abline(h=beh[k],lty=k,col=k,lwd=1.5)
#       }
#       
#       legend("topright",names,lty=1:K,col=1:K)
#       
#       vs = lapply(res,function(x){x$vn.betagammalambda})
#       vs = unlist(vs)
#       sigma.vs = lapply(res,function(x){x$sigam.vn.})
#       sigma.vs = unlist(sigma.vs)
#       zse.v = 1.96*sigma.vs/sqrt(n)
#       
#       
#       plot(lambdas,vs,ylab=TeX('$V_{n}(train)$'),xlab=TeX("$\\lambda$"),type='l',
#            ylim=c(min(c(vs-zse.v,vs)),max(c(vs+zse.v,vs))))
#       
#       lines(lambdas,vs-zse.v,lty=2)
#       polygon(c(rev(lambdas),lambdas),c(rev(vs-zse.v),
#                                         vs+zse.v),
#               col=adjustcolor(1,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
#               lty=0)
#       lower.bound = vs-zse.v
#       
#       star.ix = which(lower.bound==max(lower.bound))
#       points(lambdas[star.ix],lower.bound[star.ix],pch=8)
#       
#       star.ix.v = which(vs==max(vs))
#       points(lambdas[star.ix.v],vs[star.ix.v],pch=8)
#       
#       if(!is.null(or.shell.var)){
#         vars.n.z = unlist.reshape(lapply(res.var, function(x){(x$vn.betagammalambda)}))
#         sigma.vs.n = sqrt(vars.n.z)
#         zse.v.n = 1.96*sigma.vs.n/sqrt(M)
#         lines(lambdas,vs-zse.v.n,lty=3,col="green")
#       }
#       
#       
#       plot(lambdas,sigma.vs,ylab=TeX('$\\sigma(V_{n})(train)$'),xlab=TeX("$\\lambda$"),type='l',
#       )  
#       
#       #scale.plot(lambdas,sigma.vs,(vs[1]-vs),"train") #function(ls,x,y,label){
#       
#       ## test
#       
#       vs = lapply(res,function(x){x$vn.betagammalambda.test})
#       vs = unlist(vs)
#       sigma.vs = lapply(res,function(x){x$sigam.vn.test})
#       sigma.vs = unlist(sigma.vs)
#       zse.v = 1.96*sigma.vs/sqrt(n)
#       
#       plot(lambdas,vs,ylab=TeX('$V_{n}(test)$'),xlab=TeX("$\\lambda$"),type='l',
#            ylim=c(min(c(vs-zse.v,vs)),max(c(vs+zse.v,vs))))
#       
#       
#       
#       lines(lambdas,vs-zse.v,lty=2)
#       polygon(c(rev(lambdas),lambdas),c(rev(vs-zse.v),
#                                         vs+zse.v),
#               col=adjustcolor(1,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
#               lty=0)
#       lower.bound = vs-zse.v
#       star.ix = which(lower.bound==max(lower.bound))
#       points(lambdas[star.ix],lower.bound[star.ix],pch=8)
#       #text(lambdas[star.ix],lower.bound[star.ix],"V-2sigma/sqrt(n)")
#       
#       star.ix.v = which(vs==max(vs))
#       points(lambdas[star.ix.v],vs[star.ix.v],pch=8)
#       # text(lambdas[star.ix],vs[star.ix],"V",offset=.01)
#       
#       if(!is.null(or.shell.var)){
#         vars.n.z = unlist.reshape(lapply(res.var, function(x){(x$vn.betagammalambda.test)}))
#         sigma.vs.n = sqrt(vars.n.z)
#         zse.v.n = 1.96*sigma.vs.n/sqrt(n)
#         lines(lambdas,vs-zse.v.n,lty=3,col="green")
#       }
#       
#       plot(lambdas,sigma.vs,ylab=TeX('$\\sigma(V_{n}(test))$'),xlab=TeX("$\\lambda$"),type='l',
#       )  
#       
#       #scale.plot(lambdas,sigma.vs,(vs[1]-vs),"test") #function(ls,x,y,label){
#       
#     }
#   }
# }

get.cov=function(theor.vars,betans.m,beta0,n){
  K=length(beta0)
  M=length(theor.vars)
  theor.vars = unlist.reshape(theor.vars)
  theor.var.mean=apply(theor.vars,2,mean)
  
  quant=qnorm(0.975)
  mult = sqrt(theor.var.mean/n)*quant
  
  ucb = t(apply(betans.m,1,function(x){x+mult}))
  ucb = matrix(ucb,nrow=M,ncol=length(beta0))
  lcb =  t(apply(betans.m,1,function(x){x-mult}))
  lcb = matrix(lcb,nrow=M,ncol=length(beta0))
  
  
  #t(apply(lcb,1,function(x){beta0>=x}))
  #t(apply(ucb,1,function(x){beta0<=x}))
  covs = matrix(nrow=M,ncol=K)
  for (j in 1:length(beta0)){
    for (i in 1:M){
      covs[i,j]=((lcb[i,j]<=beta0[j]) & (beta0[j]<=ucb[i,j]))
    }
  }
  cov=apply(covs,2,mean)
  list(cov=cov,theor.var.mean=theor.var.mean,ci.width=mult*2)
}

get.cov.table =  function(param){

  deltas = param$deltas
  gammas = param$gammas
  ngam = length(gammas)
  b0 = param$b0
  n = param$n
  tau=param$tau
  T=param$T
  R=param$R
  M= param$M
  K = param$K
  sd.epsi=param$sd.epsi
  mean.first=param$mean.first
  sd.first=param$sd.first
  u=param$u
  plt.obj=param$plot.obj
  scale.s=param$scale.s
  mysel=param$mysel
  other.args = param$other.args
  
  t0=Sys.time()

  ii = 1
  mts = list()
  #for (d in 1:length(deltas)){
  d = 1
  delta=deltas[d]
  inner.list = list()
  for (g in 1:length(gammas)){
    # for each gamma, we want identical simulated datasets
    set.seed(1)
    gamma = gammas[g]
    print("need to make sim.db take other arguments")
    #set.seed(seed)
    exps = list()
    for (j in 1:M){
      if (j%%20 == 0){
        print(c("MC",j))
      }
      #para=readRDS(paste0("param/",ii,paste0(exp.tag,"exp.param")))
      #eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau,R)
      eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,sd.epsi=sd.epsi,
                 mean.first=mean.first,sd.first=sd.first,seed=j)
      
      # note that gamma_{n,gamma,lambda} cov is not really 
      # assessed, so lambda n doesn't matter
      lambda=lambdan(n,delta)
      exps[[j]]=process.eps(eps=eps,b0=b0,gamma=gamma,
                            u=u,delta=delta,lambda=lambda,
                            plt.obj=plt.obj,scale.s=scale.s,
                            sel=mysel,split=0,other.args=other.args)
      
      ii=ii+1
    }
    print("if gamma =1")
    print(c("n",n))
    print(c("mysel",mysel))
    tilde.betans = unlist.reshape(lapply(exps,function(x){x$tilde.betan}))
    print("beta_{n,gamma,0} var")
    print(apply(tilde.betans,2,var))
    betans = unlist.reshape(lapply(exps,function(x){x$betan}))
    print("beta_{n} var")
    print(apply(betans,2,var))
    
  
    #saveRDS(exps,paste0(exp.tag,"ms",ms,"exps"))
    
    #asses.cov = function(exps){
    #n = exps[[1]]$n.train
    betans = lapply(exps,function(x){x$betan})
    betanlambdas = lapply(exps,function(x){x$betanlambda})
    bns = lapply(exps,function(x){x$bn})
    #bns.optim = lapply(exps,function(x){x$bn.optim})
    tilde.betans = lapply(exps,function(x){x$tilde.betan})
    
    betans.m = unlist.reshape(betans)

    betanlambdas.m = unlist.reshape(betanlambdas)
    bns.m = unlist.reshape(bns)
    
    bns.optim = lapply(exps,function(x){x$bn.optim})
    #bns.optim.m = unlist.reshape(bns.optim)
    
    rnbetans.m = sqrt(n)*betans.m
    rnbetanlambdas.m = sqrt(n)*betanlambdas.m
    rnbns.m = sqrt(n)*bns.m
    
    varns=apply(rnbetans.m,2,var)

    varns.lambda = apply(rnbetanlambdas.m,2,var)
    varns.b = apply(rnbns.m,2,var)
    betans.mean=apply(betans.m,2,mean)
    betanslambdas.mean = apply(betanlambdas.m,2,mean)
    bns.mean = apply(bns.m,2,mean)
    
    tilde.betans.m = unlist.reshape(tilde.betans)
    beta0=apply(tilde.betans.m,2,mean)
    
    theor.vars = lapply(lapply(exps, function(x){x$var}),function(z){diag(z)})
    

    theor.vars.plugintrue = lapply(lapply(exps, function(x){x$var.plugintrue}),function(z){diag(z)})

    theor.vars.lambda = lapply(lapply(exps, function(x){x$var.betanlambda}),function(z){diag(z)})
    theor.vars.b = lapply(lapply(exps,function(x){x$var.b}),function(z){diag(z)})
    

    plot.data=1
    if (plot.data){
    main.title = paste(expression(gamma),"=",gamma,"; n=",n,"; T=",T)
    png(paste0(main.title,"cov.plot.png"),res=100)
    par(mfrow=c(2,1))
    hist(betans.m,main=main.title)
    hist(unlist(theor.vars),main=paste0("target sigma^2_0=",varns))
    #hist(theor.vars.b)
    dev.off()
    }
    
    betan.cov=get.cov(theor.vars,betans.m,beta0,n)
    betan.cov.plugintrue = get.cov(theor.vars.plugintrue,betans.m,beta0,n)
    betan.lambda.cov = get.cov(theor.vars.lambda,betanlambdas.m,beta0,n)
    betan.lambda.bcov =  get.cov(theor.vars.b,betanlambdas.m,b0,n)
    bn.bcov = get.cov(theor.vars.b,bns.m,c(b0),n)
    
    mt=rbind(
      beta0,
      betans.mean,
      (betans.mean-beta0),
      betanslambdas.mean,
      bns.mean,
      c(b0),
      sqrt(betan.cov$theor.var.mean),
      sqrt(betan.cov.plugintrue$theor.var.mean),
      sqrt(varns),
      sqrt(betan.lambda.cov$theor.var.mean),
      sqrt(varns.lambda),
      sqrt(bn.bcov$theor.var.mean),
      sqrt(varns.b),
      betan.cov$cov,
      betan.cov$ci.width,
      NA,#c(betan.lambda.bcov$cov[1],betan.lambda.cov$cov[2]),
      bn.bcov$cov
      #sqrt(betan.cov$theor.var.mean/n)*qnorm(0.975)
    )
    # bias
    
    mt = matrix(mt[,1:K],nrow=dim(mt)[1],ncol=K)
   
    colnames(mt)= paste0("$S_",1:K,"$")
    rownames(mt)=c("$\\bar{\\beta}_{0,n,\\gamma}$",
                   "$\\bar{\\beta}_{n,\\gamma}$",
                   "Bias",
                   "$\\bar{\\beta}_{n,\\gamma,\\lambda}$",
                   "$\\bar{b}_n$",
                   "$b_0$",
                   "$\\bar{\\sigma}_n(\\beta_{n,\\gamma})$",
                   "$\\bar{\\sigma}_n(\\beta_{0,n,\\gamma})$",
                   "$\\sigma_{0,n}(\\beta_{n,\\gamma})$",    
                   "$\\bar{\\sigma}(\\beta_{n,\\gamma,\\lambda})$",
                   "$\\sigma_n(\\beta_{n,\\gamma,\\lambda})$",
                   "$\\bar{\\sigma}(b_n) $",
                   "$\\sigma_n(b_n)$",               
                   "$\\beta_{n,\\gamma}$ coverage",
                   "Length CI",
                   #"Cov.",
                   "$\\beta_{n,\\gamma,\\lambda}$ coverage",
                   "$b_n$ coverage")
    
    #remove adaptive lasso
    #lo = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE)
    
    #ftt2=t(ftt[lo,])
    print("add length CI")
    lo  = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)
    cbind(rownames(mt),lo)
    mt = mt[lo,,drop=FALSE]

    mt = formatC(mt,digits=2,format="f")
    #mt=as.data.frame(mt)
    
    #if (K==1){
    #if(class(mt)=="character"){
    #  #mt = as.data.frame(mt)
    #  mt = as.matrix(mt)
    #}
    #}
    print(xtable(mt, caption = paste0("$n=$",n,", $T=$",T,", $K=$",K, ", $M=$",M,", 
$\\gamma$=",gamma, ", $\\delta$=",delta)),
          type="latex",sanitize.text.function = function(x){x})
    
    inner.list[[g]] = list(mt=mt,par=list(gamma=gamma,delta=delta))
  }
  
  mts[[d]] = inner.list
  mmts = list(mts=mts,outer.par = list(n=n,T=T,K=K,M=M))
  #saveRDS(mmts,paste0(exp.tag,"fttt.table"))
  #}
  
  print("Need to deal with SCALING")
  print("can scaling just be with na.rm??? n in other places...")
  print("NOTE R AND EST.VO change with K")
  
  at=mmts$mts[[1]][[1]]$mt
  ft = matrix(NA,nrow=dim(at)[1]*1,ncol=dim(at)[2]*length(gammas))
  rseq=seq(1,dim(ft)[1],by=dim(at)[1])
  cseq=seq(1,dim(ft)[2],by=dim(at)[2])
  
  d= 1
  #for (d in 1:length(deltas)){
  #browser()
  

  for (g in 1:length(gammas)){
    ft[rseq[d]:(rseq[d]+dim(at)[1]-1),cseq[g]:(cseq[g]+dim(at)[2]-1)] = mmts$mts[[d]][[g]]$mt
  }
  #}
  ft = as.data.frame(ft)
  
  delta.row = round(rep(deltas,each=dim(at)[1]),4)
  delta.row = formatC(delta.row,digits=2,format="f")
  #gamma.col=round(rep(gammas,each=dim(at)[2]),4)
  gamma.col = formatC(rep(gammas,each=dim(at)[2]), format = "e", digits = 2)
  #wc = cbind(deltas=delta.row,ft)
  #wc
  ftt = rbind(c(gamma.col),ft)
  #ftt[1,"rownames"] = "$\\gamma$"
  ftt
  
  # removes variables that weren't selected
  ftt2=ftt[,as.logical(c(rep(mysel,ngam)))]
  #colnames(ftt) = c("",rep(colnames(ftt)[2],dim(ftt)[2]-1))
  #rownames(ftt) = c("$\\gamma$",rownames(ftt))
  
  ftt2 = as.matrix(ftt2)

  #colnames(ftt2) = c(rep(colnames(at),length(gammas)))[as.logical(c(rep(mysel,ngam)))]
  ftt2 = rbind(c(rep(colnames(at),length(gammas)))[as.logical(c(rep(mysel,ngam)))],ftt2)
  #colnamesft = rep(colnames(at),length(gammas))

  rownames(ftt2) = c("Covariate","$\\gamma$",rep(rownames(at),1))

  ftt2 = ftt2[c("$\\gamma$", "Covariate",                           
                "$\\bar{\\beta}_{0,n,\\gamma}$",         
              "$\\bar{\\beta}_{n,\\gamma}$",            
              "Bias",
              "$\\bar{\\sigma}_n(\\beta_{n,\\gamma})$",
              "$\\bar{\\sigma}_n(\\beta_{0,n,\\gamma})$",
              "$\\sigma_{0,n}(\\beta_{n,\\gamma})$",  
              #"$\\bar{\\sigma}(\\beta_{n,\\gamma})$",  
              #"$\\bar{\\sigma}(\\beta_{0,n,\\gamma})$", 
              #"$\\sigma_n(\\beta_{n,\\gamma})$",        
              "$\\beta_{n,\\gamma}$ coverage",         
              "Length CI"),,drop=FALSE]
  #ftt2 = ftt2[c(2,1,3,4,5,6,7,8,9),]
  #saveRDS(list(tab=ftt2,ms=ms),paste0(sel.tab.dir,"/",ms))
   
  print(xtable(ftt2, caption = paste0("$n=$",n,", $T=$",T,", $K=$",K, ", $M=$",M,", Active=",
                                      paste0(mysel,sep=" ",collapse=""))),
        type="latex",sanitize.text.function = function(x){x},
        include.rownames=TRUE, include.colnames=FALSE,
        digits=2)
  print(c("sel",mysel,", Time (minutes)",(Sys.time()-t0)/60))
  
  ftt2
}

#make.post.select.table = function(coef.table){
#   beta.and.ci=paste0(coef.table$est," (",coef.table$LCI,", ",coef.table$UCI,")")
#   b.and.ci = paste0(formatC(ests.b,format="f",digits=rdigits)," (",formatC(ests.b-cbs.b,format='f',digits=rdigits),
#                     ", ",formatC(ests.b+cbs.b,format="f",digits=rdigits),")")
#   newcoef.table=cbind(beta.and.ci,b.and.ci,select.mask)
#   colnames(newcoef.table)=c("$\\beta_{n,\\gamma}$","$b_n$","Active")
#   rownames(newcoef.table)= outer.res$names
#   print(newcoef.table)
#   newcoef.table = as.data.frame(newcoef.table)
#   
#   #print(xtable(newcoef.table[newcoef.table[,'Active']==1,1:2], caption = paste0("$n=$",n.sel,", $\\gamma=$",
#   newcoef.table.set = newcoef.table
#   newcoef.table.set[newcoef.table.set[,'Active']==0,1] = "set to $b_n$"
#   
#   print(xtable(newcoef.table.set[,1:2], caption = paste0("Held out dataset post-selection inference. $n=$",n.sel,", $\\gamma=$",
#                                                          #round(gammas[selection.ix],4)
#                                                          #formatC(gammas[selection.ix], format = "e", digits = 2), ", $\\delta=$",deltas[d])),
#                                                          formatC(gammas[selection.ix], format = "e", digits = 2))),
#         type="latex",sanitize.text.function = function(x){x},include.rownames=TRUE,digits=4)
# }

make.gamma.tables = function(outer.res.sel,rdigits=2,d,n.sel,select.mask,nolam,M=NULL){

  gammas = outer.res.sel$gammas
  gamma.tables = list()
  d=1 #gamma tables dnd delta. old version had deltas there
  for (selection.ix in 1:length(gammas)){
    cbs=sqrt(diag(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$var.betanlambda))/sqrt(n.sel)*qnorm(0.975)
    cbs.b = sqrt(diag(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$var.b))/sqrt(n.sel)*qnorm(0.975)
    ests= c(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$betan)
    ests.b = c(outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$bn)
    
    coef.table = cbind(beta=ests,lcb=(ests-cbs),
                       ucb=(ests+cbs),sel=select.mask,
                       bn=outer.res.sel$outer.res[[d]][[selection.ix]][[1]]$bn)
    rownames(coef.table)= outer.res.sel$names #paste0("$S_",1:K,"$")
    colnames(coef.table)=c("est", "LCI","UCI","Active","$b_n$")
    print("coef.table")
    print(coef.table)
    if (sum(coef.table[,"Active"])==0){
      print("No active coefficients, going to break coef table I guess bc no selection")
      print(c("gamma ix=",selection.ix))
  
    }
    coef.table=formatC(coef.table,format='f',digits=rdigits)
    coef.table=as.data.frame(coef.table)
    gamma.tables[[selection.ix]] = coef.table
    beta.and.ci=paste0(coef.table$est," (",coef.table$LCI,", ",coef.table$UCI,")")
    b.and.ci = paste0(formatC(ests.b,format="f",digits=rdigits)," (",formatC(ests.b-cbs.b,format='f',digits=rdigits),
                      ", ",formatC(ests.b+cbs.b,format="f",digits=rdigits),")")
    newcoef.table=cbind(beta.and.ci,b.and.ci,select.mask)
    colnames(newcoef.table)=c("$\\beta_{n,\\gamma}$","$b_n$","Active")
    rownames(newcoef.table)= outer.res.sel$names
    print(newcoef.table)
    newcoef.table = as.data.frame(newcoef.table)
    
    #print(xtable(newcoef.table[newcoef.table[,'Active']==1,1:2], caption = paste0("$n=$",n.sel,", $\\gamma=$",
    newcoef.table.set = newcoef.table
    newcoef.table.set[newcoef.table.set[,'Active']==0,1] = "set to $b_n$"
    
    print(xtable(newcoef.table.set[,1:2], caption = paste0("Held out dataset post-selection inference. $n=$",n.sel,
                                                           ", $M=$",M,", $\\gamma=$",
                                                           #round(gammas[selection.ix],4)
                                                           #formatC(gammas[selection.ix], format = "e", digits = 2), ", $\\delta=$",deltas[d])),
                                                           formatC(gammas[selection.ix], format = "e", digits = 2))),
          type="latex",sanitize.text.function = function(x){x},include.rownames=TRUE,digits=4)
    
    #make.post.select.table(coef.table)
      beta.and.ci=paste0(coef.table$est," (",coef.table$LCI,", ",coef.table$UCI,")")
      b.and.ci = paste0(formatC(ests.b,format="f",digits=rdigits)," (",formatC(ests.b-cbs.b,format='f',digits=rdigits),
                        ", ",formatC(ests.b+cbs.b,format="f",digits=rdigits),")")
      newcoef.table=cbind(beta.and.ci,b.and.ci,select.mask)
      colnames(newcoef.table)=c("$\\beta_{n,\\gamma}$","$b_n$","Active")
      rownames(newcoef.table)= outer.res.sel$names
      print(newcoef.table)
      newcoef.table = as.data.frame(newcoef.table)

      #print(xtable(newcoef.table[newcoef.table[,'Active']==1,1:2], caption = paste0("$n=$",n.sel,", $\\gamma=$",
      newcoef.table.set = newcoef.table
      newcoef.table.set[newcoef.table.set[,'Active']==0,1] = "set to $b_n$"

      print(xtable(newcoef.table.set[,1:2], caption = paste0("Held out dataset post-selection inference. $n=$",n.sel,", $\\gamma=$",
                                                             #round(gammas[selection.ix],4)
                                                             #formatC(gammas[selection.ix], format = "e", digits = 2), ", $\\delta=$",deltas[d])),
                                                             formatC(gammas[selection.ix], format = "e", digits = 2),
                                                            ", No $\\lambda$ satisfied=",nolam)
                   ),
            type="latex",sanitize.text.function = function(x){x},include.rownames=TRUE,digits=4)
  }
}

run.one.sel = function(param,m){
  exp.tag = param$exp.tag
  
  resfile = paste0(exp.tag,"/lambplots.outerress/",m)
  
  deltas = param$deltas
  gammas = param$gammas
  lambdas = param$lambdas
  ngam = length(gammas)
  ndel = length(deltas)
  nlam = length(lambdas)
  b0 = param$b0
  n = param$n
  tau=param$tau
  T=param$T
  R=param$R
  M= param$M
  K = param$K
  sd.epsi=param$sd.epsi
  mean.first=param$mean.first
  sd.first=param$sd.first
  u=param$u
  plt.obj=param$plt.obj
  scale.s=param$scale.s
  mysel=param$mysel

  res.dir.sel = param$res.dir.sel
  rdigits = param$rdigits
  names =  param$names
  other.args = param$other.args
  
  gamma.select.ix = other.args$gamma.select.ix
  lambda.select.ix = param$lambda.select.ix
  #delta.select.hardcode=1,gamma.select.ix=1
  
  eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,
             sd.epsi=sd.epsi,mean.first=mean.first,sd.first=sd.first,
             seed=m)
  plt.obj=FALSE
  if (plt.obj){
    plot.objs(eps,gammas,u,b0)
    browser()

    
  }
  


  n = length(eps)
  n.train = as.integer(length(eps)/2)
  n.test = n-n.train
  #n.test = outer.res$outer.res[[1]][[1]][[1]]$n.test
  #n.train =  outer.res$outer.res[[1]][[1]][[1]]$n.train
  #n=outer.res$outer.res[[1]][[1]][[1]]$n.test
  train.eps = eps[1:n.train]
  test.eps=eps[(n.train+1):n]
  
  
  outer.res = lambda.path(train.eps,b=b0,gammas=gammas,lambdas=lambdas,
                          names=names,deltas=deltas,
                          file="mc",split=1,scale.s=0,sel=rep(1,K),
                          lambda.select.ix=lambda.select.ix,
                          gamma.sel=gamma.select.ix,other.args=other.args)
  
  
  
  ixs=plot.lambda.paths(outer.res)
  ##
  gamma.select.ix = other.args$gamma.select.ix
  delta.hardcode = other.args$delta.select.hardcode#2
  #browser()
  diff.gt.z = ixs[[delta.hardcode]][[ngam]]$roll.diff.p.gt.z

  #if (lambda.select.ix){
  #  select.mask = lambda.select.ix
  #}else{
  select.mask = ixs[[delta.hardcode]][[gamma.select.ix]]$selected.mask
  #}
  print(ixs)
  nolam = ixs[[delta.hardcode]][[gamma.select.ix]]$no.lambda.satisfies
  print(c("No lambda satisfies:",nolam))
  print("select mask")
  print(select.mask)
  # thought about making a table for each possible one
  #possible.selects = get.possible.selects(K,1)
  #outer.res.sels = list()
  #for (sm in 1:length(possible.selects)){
 
    outer.res.sel = lambda.path(test.eps,b0=b0,gammas=gammas,lambdas=c(0),
                                names=names,
                                #names=c(TeX('$S_1$'),TeX('$S_2$')),
                                deltas,file="mc",
                                sel=select.mask,split=0,scale.s=scale.s,
                                lambda.select.ix=lambda.select.ix,
                                gamma.sel=gamma.select.ix, other.args=other.args)
    
   # outer.res.sels[[sm]] = outer.res.sel
    # note that n comes from now sel, which should be just half data
    n.sel=outer.res.sel$n
    
    print(c("correct n (check)?",n.sel==n.test))
    # we only take half the data
    
    # selection ix corresponds to which gamma we select for
    #selection.ix = gamma.select.ix
    
    # it's outer.res[[delta]][[gamma]][[lambda]]. lambda is 1, which is 0, since we are doing this
    # after selection.
    
    #for (d in 1:length(deltas)){
    # does not depend on delta
    
    d=1
   
    gamma.tables = make.gamma.tables(outer.res.sel,rdigits,d,n.sel,select.mask,nolam)
    # }

  
  #}
  oplist = list(outer.res=outer.res, 
                #outer.res.sels=outer.res.sels,
                outer.res.sel=outer.res.sel,
                gamma.tables=gamma.tables,
                select.mask=select.mask,
                diff.gt.z=diff.gt.z)
  saveRDS(oplist,paste0(res.dir.sel,"/",m))
  
}

list.outer.res.to.summ = function(ors,mean=1){
  
  outer.res = ors[[1]]#$outer.res
  or.shell.av = outer.res
  or.shell.var = outer.res
  gammas = outer.res$gammas
  deltas = outer.res$deltas
  lambdas = outer.res$lambdas
  # outer.res$outer.res is [[delta]][[gamma]][[lambda]]
  
  # all lambdas. so i think we would be doing unlist(ors[[1]]$outer.res 
  # to average over all delta gamma lambda instead of just lamdba
  print(c("outer.res",outer.res))
  number.lam.gam.del.items = length(unlist(outer.res$outer.res))
  # for each MC iteration, make a column.  
  # The rows will be the number of items in one of the 
  # results lists (indexed by a gamma,delta, and lambda)
  M =length(ors)
  r = matrix(nrow=number.lam.gam.del.items,ncol=M)
  for (m in 1:M){
    r[,m] = unlist(ors[[m]]$outer.res)
  }
  
  # then run stats. note that even taking averages 
  # of like n.train, which is just n.train
  #for (mean in c(0,1)){
  mean = 1
  if (mean){
    rav = apply(r,1,mean)
  }else{
    rav = apply(r,1,median)
  }
  rvar = apply(r,1,var)
  
  # one lambda
  # just get the names of one element of outer.res (one process.eps output)
  m.names=names(outer.res[[1]][[1]][[1]][[1]])
  k=1
  for(z in 1:length(deltas)){
    for (j in 1:length(gammas)){
      for (i in 1:length(lambdas)){
        for (na in m.names){
          
          # get the dimension of this named item (just use delta=gamma=lambda=1, since 
          # dimension is the same for all delta, gamma, lambda)
          # the last [[1]] is just to get the item itself
          m.dim = dim(as.matrix(or.shell.av$outer.res[[1]][[1]][[1]][na][[1]]))
          # cast as matrix so can get dim.  In the plot.lambda paths, we will uncast
          # to scalar again for some, such as sample size n
          if(is.null(m.dim)){m.dim=c(1,1)}
          # now because we unlisted this item, it should be rav dimensions 
          # k:k+(m.dim[1]*m.dim[2]-1), minus one is there because we start with k=1 not k=0
          av.res = matrix(rav[k:(k+(m.dim[1]*m.dim[2]-1))],
                          nrow=m.dim[1],ncol=m.dim[2])
          va.res = matrix(rvar[k:(k+(m.dim[1]*m.dim[2]-1))],
                          nrow=m.dim[1],ncol=m.dim[2])
          or.shell.av$outer.res[[z]][[j]][[i]][na][[1]] = av.res
          or.shell.var$outer.res[[z]][[j]][[i]][na][[1]] = va.res
          k=k+m.dim[1]*m.dim[2]
        }
      }
    }
  }
  ### ADDED 
  or.shell.av$or.shell.var = or.shell.var
  
  ### ADDED
  or.shell.av
}

outer.res.to.arr = function(ors.av,metric,trans=function(x){x}){
  deltas = ors.av$deltas
  ndel = length(deltas)
  gammas = ors.av$gammas
  ngam = length(gammas)
  lambdas = ors.av$lambdas
  nlam = length(lambdas)
  sear.matr = array(NA,c(ndel,ngam,nlam),dimnames=list(deltas=deltas,gammas=gammas,lambdas=lambdas))
  for (d in 1:ndel){
    for (g in 1:ngam){
      for (l in 1:nlam){
        sear.matr[d,g,l] = trans(ors.av$outer.res[[d]][[g]][[l]][[metric]])
      }
    }
  }
  sear.matr
}

analyze.sel.res = function(res.dir,sel.res.dir,plot.dir,rdigits=3){
  
  print(sel.res.dir)
  print(plot.dir)
  fs = list.files(path=sel.res.dir,full.names=TRUE)
  print(fs)
  fs = mixedsort(fs)
  ors.and.tbs = lapply(fs,readRDS)
  
  masks = lapply(ors.and.tbs,function(x){x$select.mask})
  diff.gt.zs =  lapply(ors.and.tbs,function(x){x$diff.gt.z})
  unique.masks = unique(masks)#masks[!duplicated(lapply(masks, sort))]

  #browser()
  print("Are all selected masks the same?")

  mm=unlist.reshape(masks)
  print(mm)
  print(apply(mm,2,var)==0)
  saveRDS(mm,paste0(res.dir,'/masks.matrix'))
  
  # dm=unlist.reshape(masks)
  # print(dm)
  # print(apply(dm,2,var)==0)
  # saveRDS(dm,paste0(res.dir,'/diff.gt.z.matrix'))
  
  ms = matrix(apply(mm,2,mean),nrow=1)
  m.names=ors.and.tbs[[1]]$outer.res$names 
  colnames(ms) = m.names
  rownames(ms) = "$\\hat{P}(selected)$"
  print(xtable(ms),type="latex",sanitize.text.function = function(x){x})
  
  #ds = matrix(apply(dm,2,mean),nrow=1)
  #d.names=gammas
  #colnames(ds) = d.names
  #rownames(ds) = "$\\hat{P}(|\\pi_{\\beta_{n,\\gamma}}-\\pi_{b_{n}}|>\\zeta)$"
  #print(xtable(ds),type="latex",sanitize.text.function = function(x){x})
  
  
  #ms = apply(mm,2,mean)
  #print(xtable(ms))
  
  ors = lapply(ors.and.tbs,function(x){x$outer.res})
  ors.sel =  lapply(ors.and.tbs,function(x){x$outer.res.sel})
  g.tbs =  lapply(ors.and.tbs,function(x){x$gamma.tables})
  
  M = length(ors)

  # get info to make plots
  outer.res = ors[[1]]
  
  gamma.select.ix=outer.res$other.args$gamma.select.ix
  K=dim(outer.res$outer.res[[1]][[1]][[1]]$betan)[1]
  gammas = outer.res$gammas
  ngam = length(gammas)
  deltas = outer.res$deltas
  ndel = length(deltas)
  lambdas = outer.res$lambdas
  nlam = length(lambdas)
  n = outer.res$n
  T= outer.res$T
  use.mean = 1

  ors.av = list.outer.res.to.summ(ors,mean=use.mean)
  
  
  res = ngam*res.multiplier
  nplots=3
  tag = "noTag"


  #plotfile = paste0(sel.res.dir,"/",tag,"mean",mean,
  # "_seed=",sd,"n=",n,"M=",M,"K=",K,"T=",T,"_MC.reps.png")
  plotfile = paste0(plot.dir,"/",tag,"use.mean",use.mean,
                    "n=",n,"M=",M,"K=",K,"T=",T,"_MC.reps.png")
  print(plotfile)
  #png(plotfile,height=600*ngam*nplots,width=850*ndel,res=res)
  png(plotfile,height=heightchunk*ngam*nplots,width=widthchunk*ndel,res=res)
  #par(mfcol=c(ndel*nplots,ngam)) #vertical in gamma
  

  ml=layout(matrix(c(1:(ndel*nplots*ngam)),ncol=ndel), 
            widths=rep(myw,nplots*ngam*ndel), 
            heights=rep(myheights,ngam), TRUE)
  
  #par(mfcol=c(ngam*nplots,ndel)) #horiz in gamma
  # : bottom, left, top, and right.
  par(mar=c(4.1,5.1,2,.8))

  ixs=plot.lambda.paths(ors.av,gammaixleg=1,deltaixleg=3)
  mar=get.stats(ors.av)
  saveRDS(mar,paste0(res.dir,'/stats'))

  dev.off()
  saveRDS(ors.av,paste0(res.dir,'/ors.av'))


  mat.mean.pi = outer.res.to.arr(ors.av,'rollout.betangammalambda.mean')
  mat.mean.diff = outer.res.to.arr(ors.av,'roll.diff.betangammalambda.bn.ab')
  saveRDS(mat.mean.pi,paste0(res.dir,'/mat.mean.pi'))
  saveRDS(mat.mean.diff,paste0(res.dir,'/mat.mean.diff'))
  print(mat.mean.pi)
  print(mat.mean.diff)
  
  # use the selection from the average of lambda plots
  #mysel = ixs[[gamma.select.ix]]$selected.mask
  #nolam =  ixs[[gamma.select.ix]]$no.lambda.satisfies
  nolam = NULL
  for (mysel in unique.masks){
    
  mask.indices = unlist(lapply(masks,function(x){all(x==mysel)}))

  # we only select the ones with the same mask indices
  ors.sel.av = list.outer.res.to.summ(ors.sel[mask.indices])
 
  # masks 
  n = ors.sel.av$n # or is it n.test? yes, it's 1/2 n


  make.gamma.tables(ors.sel.av,rdigits,d,n,mysel,nolam,M=sum(mask.indices)) 
  }
  saveRDS(unique.masks,paste0(res.dir,'/unique.masks'))
  unique.masks
  
}

get.stats = function(outer.res){

  other.args = outer.res$other.args
  lambdas=outer.res$lambdas
  deltas=outer.res$deltas
  gammas=outer.res$gammas
  nlam=length(lambdas)
  ndel=length(deltas)
  ngam=length(gammas)
  items = c("roll.diff.p.gt.z",
            "rollout.betangammalambda.mean",
            "mean.lp.rollout.betanlambda.sel",
            "mean.lp.rollout.betanlambda.nosel",
            "nsel")
  mar = array(NA,c(ndel,ngam,nlam,length(items)))
  for (d in 1:ndel){
    for (g in 1:ngam){
      for (l in 1:nlam){
        mitems = rep(NA,length(items))
        for (i in 1:length(items)){
          mitems[i] = outer.res$outer.res[[d]][[g]][[l]][[items[i]]]
        }
        mar[d,g,l,] = mitems
      }
    }
  }
list(mar=mar,items=items)
  
}


get.ratios = function(outer.res){
  
  other.args = outer.res$other.args
  lambdas=outer.res$lambdas
  deltas=outer.res$deltas
  gammas=outer.res$gammas
  nlam=length(lambdas)
  ndel=length(deltas)
  ngam=length(gammas)

  length.rat = length(outer.res$outer.res[[1]][[1]][[1]][["ratios"]])
  mar = array(NA,c(ndel,ngam,nlam,length.rat))
  for (d in 1:ndel){
    for (g in 1:ngam){
      for (l in 1:nlam){
        mar[d,g,l,] = outer.res$outer.res[[d]][[g]][[l]][["ratios"]]
      }
    }
  }
  
  mar
  
}

