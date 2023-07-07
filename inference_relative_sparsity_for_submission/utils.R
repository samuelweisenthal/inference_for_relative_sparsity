source('center.scale.R')
# general helpers.
# some not used for experiments in main text


slog = function(x){log(x+1e-50)}

#expit = function(s){
#  r=exp(s)/(1+exp(s))
#  as.vector(r)
#}

# MC = function(M,b0,gamma,u,delta,init.state.mean,
#               init.state.var,trans.var,T,minvar=FALSE,seed=1,tau=tau){
#   #set.seed(seed)
#   exps = list()
#   for (j in 1:M){
#     if (j%%20 == 0) {
#       print(c("MC",j))
#     }
#     
#     eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau)
#     #eps = sim.db(b0 = b0,N=n,tau=0.1,T = T,plot. = FALSE)
#     #lines(beta.grid,v0s,lty=1)
#     #lines(beta.grid,vns,lty=2,lwd=4)
#     #plot.ob(beta.grid,mint,max,plt=FALSE)
#     lambda=lambdan(n,delta)
#     exps[[j]]=process.eps(eps=eps,b0=b0,gamma=gamma,
#                           u=u,delta=delta,lambda=lambda)
#   }
#   exps
# }

# process.eps = function(eps,b0,gamma,u,delta,lambda){
#   # too big. split up into pieces that process a dataset
#   # and pieces necessary for MC. Maybe put MC pieces into MC function
#   #so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K))
#   K = dim(eps[[1]]$Ss[[1]])[1]
#   T=length(eps[[1]]$Ss)
#   so = center.scale.s(eps)
#   aa = unlist(lapply(eps,function(x){x$As}))
#   sm = so$sm.cs
#   eps=so$e
#   
#   #a=d$a
#   #s=d$s
#   bm=glm(aa~sm-1,family='binomial')
#   b=bm$coefficients
#   dont.est.beh=0
#   if (dont.est.beh){
#     b=b0
#   }
#   dont.est.beh=0 
#   #method='Nelder-Mead' # maybe write gradient descent using l1 approx
#   method='BFGS'
#   maxit=1e6
# 
#   
#   #start at behavioral to minimize positivity viol at beginning
#   optim.res = optim(par=c(b),mn,
#                     control=list(fnscale=-1,maxit=maxit),
#                     method=method,
#                     eps=eps,b=b,gamma=gamma,u=u)
# 
#   #b.optim = optim(par=c(b),l,
#   #                control=list(fnscale=-1,maxit=maxit),
#   #                method=method,
#   #                eps=eps)
#   bn.optim = NULL#b.optim$par
#   #todo: check if optim w NT obs matches optim w N length T obs (check theory)
#   #b=bn.optim
#   betan=matrix(optim.res$par)
#   vn.betagamma = vn()
#   optim.res.lamb = optim(par=c(b),wn,
#                          control=list(fnscale=-1,maxit=maxit),
#                          method=method,
#                          b=b,eps=eps,gamma=gamma,u=u,betahat=betan,
#                          bhat=b,delta=delta,lambda=lambda)
#   betanlambda=matrix(optim.res.lamb$par)
#   
#   tilde.optim.res = optim(par=c(b),est.v0.lambda,
#                           control=list(fnscale=-1,maxit=maxit),
#                           method=method,
#                           eps=eps,gamma=gamma,u=u,b=b0)
#   tilde.betan=matrix(tilde.optim.res$par)
#   
#   var = get.var(eps,betan,b,gamma,u=u,dont.est.beh=dont.est.beh)$VAR
#   var.lambda = get.var(eps,betanlambda,b,gamma,u=u,dont.est.beh=dont.est.beh)$VAR
#   var.b = get.var(eps,betanlambda,b,gamma,u=u,dont.est.beh=dont.est.beh)$VAR.B
#   
#   
#   list(var=var,var.lambda=var.lambda,
#        var.b=var.b,
#        betan=betan,betanlambda=betanlambda,
#        b=b,bn.optim=bn.optim,
#        tilde.betan=tilde.betan)
# }

d.to.eps=function(d){
  
  n=dim(d$s)[1]
  T = length(d$a)
  eps=list()
  for (i in 1:n){
    Ss=list()
    As=list()
    Rs=list()
    for (t in 1:T){
      Ss[[t]]=d$s[i,t]
      As[[t]]=d$a[i,t]
      Rs[[t]]=d$r[i,t]
    }
    eps[[i]]=list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
}

get.R.by.arg = function(R,s,a,sp){
  # sometimes easier to have 2 arg reward
  # bc then don't lose last state action pair
  # this accomodates 2 and 3 arg  
  n.args.R = length(formals(R))
  if (n.args.R==2){
    Rt = R(s,a)
  }else{
    Rt = R(s,a,sp)
  }
  Rt
}

get.Rt = function(ep,R,t){
	
  # if function.R, calculate R 
  # (eg will do so on standardized states)
  # if function.R=0, use calculated R
  # from simulation.
  # function.R=0 is better, but I had done function.R=1
  # for the experiments. I checked
  # and the results do not appear
  # to be too different between the two.
  # I think though it makes more sense to standardize x and not
  # recompute the outcome (function.R=0).
  
  function.R=0
  if (function.R){
  Rt = get.R.by.arg(R,ep$Ss[[t]],ep$As[[t]],ep$Ss[[t+1]])
  
  }else{
	# we treat R as observed 
    Rt = ep$Rs[[t]]
  }
  Rt
}

tile = function(r,dim){
  matrix(r,nrow=dim[1],ncol=dim[2],byrow=TRUE)
}



eps.log = function(x,EPS=1e-15){
  log(x+EPS)
}

argmax = function(A,b){
  mix = which.max(A%*%b)
  # really it's dist with dirac
  # could make closer to soft, same function
  # just with temp..
  list(val=A[mix,],ix=mix) 
}

softmax = function(A,b){ 
  # really should be called softargmax
  # check this
  t = 1 
  #t = 1 # actually gives near uniform
  p = exp(A%*%b/t)/sum(exp(A%*%b/t))
  list(val=t(A)%*%p,ix=p)
}

#exp(x)^100/sum(exp(x)^(100))
#exp(x)^100/sum(exp(x)^(100))*x

safe.expsumexp = function(x){
  # to prevent exp(real) = Inf, use fact that
  # softmax(x) = softmax(x+c) for all c
  x = x - max(x)
  exp(x)/sum(exp(x))
}

true.softmax = function(x){
  p = safe.expsumexp(x)
  p%*%x
}

true.softargmax = function(x){
  p = safe.expsumexp(x)
  p
}

l1norm = function(u){
  sum(abs(u))
}

real_vec_to_prob = function(real,row,phi,basis.degree){
 
  # input real vector of length N
  # output prob vector of length N+1
  #phi = function(i,j){i}
  real = c(real,0) # is this like sneakily doing optim for smaller d parameter?
  cols = 1:length(real)
  
  phis = rep(NA,length(cols))
  for (a in cols){
    phis[a] = phi(row,a,basis.degree)
  }
  real = phis*real
  true.softargmax(real)
}

prob_vec_to_real = function(prob,EPS=1e-16){
  print("Adding eps to normalizing calc. 
        check prob_vec_to_real!")
  Z = 1/(prob[length(prob)]+EPS) # this is ok, checked w vterm and lterm
  print("adding eps in log in prob_vec_to_real. 
        is this ok?")
  eps.log(prob*Z)[1:length(prob)-1] # think should add eps in case 0?
}


prob_mat_to_real = function(P){
  # why matrix() here?
  #matrix(t(apply(prob_mat,1,prob_vec_to_real)))
  if(is.null(dim(P))){browser()}
  R = matrix(nrow=dim(P)[1],ncol=(dim(P)[2]-1))
  for (row in 1:dim(P)[1]){
    R[row,] = prob_vec_to_real(P[row,])
  }
  R
}

real_dims_to_prob = function(real.dims){
  #if (length(real.dims) == 2){
  if (length(real.dims)==3){
    c(real.dims[1:2],real.dims[3]+1)
  }else if (length(real.dims)==2){
    c(real.dims[1], real.dims[2]+1)
  }else{
    print("check shape")
    browser()
  }
}

prob_dims_to_real = function(prob.dims,basis.degree="default",dim.pib.par){
    dim.pib.par
}


real_3darray_to_prob = function(beta){
  #"real_to_prob is deprecated, but still best.
  #      Have worked to update a bit so at least
  #      it is numerically safe.
  #      todo: replace with real_3darray_to_prob in utils.R.
  #      The issue is the way that slicing arrays gives vectors
  #      or something that is like a vector with 3 dim, 
  #      not matrices."
  
  Pdim = real_dims_to_prob(dim(beta))
  alpha = array(data=NA,Pdim)

  for (a in 1:Pdim[2]){
    alpha[,a,] = cbind(beta[,a,],0) #
  }
  Psasp = array(data=NA,Pdim)
  # slow, but easy to debug for now
  for (s in 1:Pdim[1]){
    for (a in 1:Pdim[2]){
      Psasp[s,a,]=safe.expsumexp(alpha[s,a,])
    }
  }
  Psasp
}


prob_3darray_to_real=function(my.prob){
  #should be abstracted into same functin as real to array
  if(is.null(dim(my.prob))){browser()}
  real.array = array(data=NA,prob_dims_to_real(dim(my.prob)))
  K = dim(my.prob)[1]
  for (k in 1:K){
    real.array[k,,] = prob_mat_to_real(my.prob[k,,])
  }
  real.array
}

flatten = function(beta){
  # dim.beta is dimension of 3d probability matrix
  # which will flatten
  # where the sums over the last index are 1
  betaf = c(beta)
  betaf
}

unflatten = function(beta.flat,dim.beta.array,basis.degree="default"){
    # does same as above anyway
  betaa = array(data=NA,dim.beta.array)
    # is it guaranteed to work?
  betaa[] = beta.flat
  betaa
}


get.p <- function(...,p){
  x <- list(...) 
    p(...)
}


get.s.stat = function(eps,stat){#dim.s,stat){
  dim.s=dim(eps[[1]]$Ss[[1]])[1]
  eplengths = unlist(lapply(eps,function(x){length(x$As)}))
  apply(
    matrix(
      unlist(
        lapply(eps, function(x){x$Ss})),
      nrow=sum(eplengths),
      ncol=dim.s,byrow=TRUE)
    ,2,stat)
}

get.a.stat = function(eps,stat){
  eplengths = unlist(lapply(eps,function(x){length(x$As)}))
  apply(
    matrix(
      unlist(
        lapply(eps, function(x){x$As})),
      nrow=sum(eplengths),
      ncol=1,byrow=TRUE)
    ,2,stat)
}

get.r.stat = function(eps,stat){
  eplengths = unlist(lapply(eps,function(x){length(x$As)}))
  apply(
    matrix(
      unlist(
        lapply(eps, function(x){x$Rs})),
      nrow=sum(eplengths),
      ncol=1,byrow=TRUE)
    ,2,stat)
}


see.v = function(eps,pi_sa,parb,gamma,R,phi,basis.degree,intercept){
  
  print("only see.v for disc a")
  minmax = 100
  gr=seq(-minmax,minmax,length.out=100)#seq(-5,5,.1)
  if (disc.a){
    gr=gr
  }else{
    gr = rbind(gr,0)
  }
  vals = c()
  vals.good.sc = c()
  vals.bad.sc = c()
  
  if(disc.a){
    lg = length(gr)
  }else{
    lg = length(gr[1,])
  }
  for (i in 1:lg){
    if (disc.a){
      par = gr[i]
    }
    else{
      par = matrix(gr[,i],nrow=dim.s+1,ncol=1)
    }
    parb2=0
    to = ope.Thomas(eps,pi_sa,par,parb2,
                    gamma,R,phi,1,intercept)
    vals[i] = to$v
    #vals.good.sc[i] = ope.Thomas(eps,pi_sa,par = ,2,
    #                     gamma,R,phi,basis.degree,intercept)
    #vals.bad.sc[i] = ope.Thomas(eps,pi_sa,par,-2,
    #                     gamma,R,phi,basis.degree,intercept)
  }
  vs=c(vals,vals.good.sc,vals.bad.sc)
  if (disc.a){
    gr=gr
  }else{
    gr = gr[1,]
  }

  plot(gr,vals,ylim=c(min(vs),max(vs)),col=1,type='l',
       ylab="hat{V}",xlab="beta")
  #real.v = function(beta){-1*(beta**2-2*beta+2)}
  #lines(gr,real.v(gr),lty=2)
  #lines(gr,vals.good.sc,col=3)
  #lines(gr,vals.bad.sc,col=2)
  
  browser()
}

runexp = function(path.to.setup,run.sim,an=TRUE){
  source(path.to.setup)
  exp.list = readRDS(paste0(save.direc,"/param"))
  
  told = Sys.time()
  exp.res = run.sim(exp.list=exp.list)
  
  print("runtime")
  rt=difftime(Sys.time(),told,units='mins')
  print(rt)
  save.direc
  if(an){
  source('analysis.R')
  }
}

plot.d.mdp = function(d){
  s=lapply(d,function(x){x$Ss})
  a = lapply(d,function(x){x$As})
  Ts =  unlist(lapply(d,function(x){length(x$As)}))
  ss=unlist.reshape.eps.ss(s)
  #aa=unlist.reshape(a)
  #if(mean){
  #  lines(1:Ts, apply(ss,2,mean))
  #}else{

  plot(1:Ts[1],ss[1,],xlim=c(0,max(Ts)),ylim=c(min(ss),max(ss)),pch=NULL,col='white')
  for (i in 1:n){
    lines(1:Ts[i],ss[i,],cex=.1,col=rgb(0,0,0,alpha=0.9))
    }
  #}
}

unlist.reshape=function(theta.stars){
  theta.stars.m = matrix(unlist(theta.stars),
                         nrow=length(theta.stars),
                         ncol=length(theta.stars[[1]]),
                         byrow=TRUE)
  theta.stars.m
}

unlist.reshape.eps.ss = function(ss){
  # Makes a n x T matrix for plotting.
  # only good for k =1
  Ts = unlist(lapply(ss,function(x){length(x)}))
  matrix(unlist(ss),nrow=length(ss),ncol=max(Ts),byrow=TRUE)
}

unlist.reshape.eps.ss.sumTsxK = function(ss){
  # makes a sum(Ts)xK matrix for scaling.
  # we unroll over time to get the variance for scaling
  # since the behavioral policy is stationary
  # and therefore sees time as stationary.
  # this is problematic when the variance of something increases over time
  # leads to positivity violation
  # hence simulations avoid this, and hopefully real data avoids this
  Ts = unlist(lapply(ss,function(x){length(x)}))
  if (is.null(ss[[1]][[1]])){
    print("Null ss")
    browser()
  }
  K=dim(ss[[1]][[1]])[1]
  matrix(unlist(ss),nrow=sum(Ts),ncol=K,byrow=TRUE)
}


