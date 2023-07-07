#has plots
library(numDeriv)
library(Deriv)
source('grad.R')
source('is.R')
source('liops.R')
source('est.v0.R')
source('sim.R')
source('utils.R')
# expit = function(s){
#   r=exp(s)/(1+exp(s))
#   as.vector(r) # bc if matrix can't scalar multiply by other matrices as easily
# }
# slog = function(x){log(x+1e-50)}
#expit = function(x){
#  # ths is ideal, but Auto diff Deriv can't diff the max

#  k  = max(0,x)#(0>x)*0+(x>=0)*x#max(0,x)
# exp(x-k)/(exp(0-k)+exp(x-k))
#}


cross.diff = list()
grads = list()
li.grads = list()
sss = list()
NC=200
for (i in 1:NC){
  set.seed(i)
gamma=0.9
Ss=c(1,2)
As=c(1,0)
Rs=c(1,2)
T=length(Ss)
u=1
#### numerical dervis (better for multistage, not exact)
b=1
beta=2
round(grad(kappa,x=beta,Ss=Ss,As=As),7)==round(dthetakappa(beta,Ss,As),7)
(round(grad(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u),7)
  ==round(dbetami(beta,b,Ss,As,Rs,gamma,u=u),7))
(round(hessian(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,u=u,Rs=Rs),7)
  ==round(ddbetami(beta=beta,b=b,Ss=Ss,As=As,Rs=Rs,u=u,gamma=gamma),7))
round(dli(Ss,As,b),7)==round(grad(li,x=b,Ss=Ss,As=As),7)
#round(ddlidbdb(1,1,1),7)==
round(ddli(Ss,As,b),3)==round(hessian(li,x=b,Ss=Ss,As=As),3)
#gnd.grad = function(b){grad(mi,x=beta,b=b,Ss=Ss,As=As,Rs=Rs,gamma=gamma,u=u)}
#nd.cross = grad(gnd.grad,x=b) # deosn't work for vector inputs.
#round(nd.cross,5)==round(dbdbetami(beta,b,Ss,As,Rs,gamma,u),5)
#######


eps=list(
  list(As=list(1,0),
       Ss=list(matrix(c(1,1)),matrix(c(4,2))),
       Rs=list(1,5)),
  list(As=list(1,1),
       Ss=list(matrix(c(1,2)),matrix(c(1,2))),
       Rs=list(2,3))
)
K=2
b0 = matrix(rep(0,K))
init.state.mean = matrix(rep(0,K))
init.state.var=diag(rep(1,K))
trans.var = diag(rep(1,K))
n=1 #actually note that I just take the first episode...


R = function(s,a,T){
#-s[2]*a
 r=-s[2]*a
#  if(dim(r)[1]>1){
#    print("R dim greater than 1")
#    browser()
 # }
 r
}

b0 = matrix(rep(0,K))
gamma=1
init.state.mean = matrix(rep(0,K))
init.state.var=.1*diag(K)
trans.var = .1*diag(K)
sd.epsi=diag(rep(1,K))
mean.first=matrix(rep(0,K))
sd.first=diag(rep(1,K))
tau=.1
#eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau=tau,R=R)

eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,
           sd.epsi=sd.epsi,mean.first=mean.first,sd.first=sd.first)
#eps = sim.db(b0=b0,N=n,tau=.1,T=10,R=R)
# just generate

Ss=eps[[1]]$Ss
As=eps[[1]]$As
Rs=eps[[1]]$Rs

sss[[i]] = unlist.reshape(Ss)[,2]
T=length(Ss)
u=1
#### numerical dervis (better for multistage, not exact)
beta= matrix(rep(1,K))#matrix(c(1,2))[1:K]
b=matrix(rep(0,K))
round(grad(kappa,x=beta,Ss=Ss,As=As),7)==round(dthetakappa(beta,Ss,As),7)
gradmi.nd = grad(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u)
gradmi = dbetami(beta,b,Ss,As,Rs,gamma,u)
round(gradmi.nd,7)==round(gradmi,7)
(round(hessian(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u),7)
  ==round(ddbetami(beta=beta,b=b,Ss=Ss,As=As,Rs=Rs,gamma=gamma,u=u),7))
grad.li = dli(Ss,As,b) 
grad.li.nd = grad(li,x=b,Ss=Ss,As=As)
li.grads[[i]]=grad.li[1]-grad.li.nd[1]
round(grad.li,7)==round(grad.li.nd,7)
#round(ddlidbdb(1,1,1),7)==
round(ddli(Ss,As,b),4)==round(hessian(li,x=b,Ss=Ss,As=As),4)
#gnd.grad = function(b){grad(mi,x=beta,b=b,Ss=Ss,As=As,Rs=Rs,gamma=gamma)}
#nd.cross = grad(gnd.grad,x=beta) # deosn't work for vector inputs.
#round(nd.cross,5)==

# the cross products are more difficult to check,
# we will have to do it one at a time
if (K>1){
cr=round(dbdbetami(beta,b,Ss,As,Rs,gamma,u),5)

beta1 = beta[1]
beta2 = beta[2]
b1 = b[1]
b2= b[2]

mi.ind.11 = function(beta1,b1){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u)
}
grad.m11 = function(b1){
  grad(mi.ind.11,x=beta1,b1=b1)
}

mi.ind.22 = function(beta2,b2){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u)
}

grad.m22 = function(b2){
  grad(mi.ind.22,x=beta2,b2=b2)
}

c22=grad(grad.m22,x=b2)

c11=grad(grad.m11,x=b1)

mi.ind.12 = function(beta1,b2){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u)
}

grad.m12 = function(b2){
  grad(mi.ind.12,x=beta1,b2=b2)
}

c12=grad(grad.m12,x=b2)

mi.ind.21 = function(beta2,b1){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u)
}

grad.m21 = function(b1){
  grad(mi.ind.21,x=beta2,b1=b1)
}

c21=grad(grad.m21,x=b1)

# this seems to depend sometimes on the seed. if not equal then close
nd.cr=matrix(c(c11,c21,c12,c22),nrow=2,ncol=2)
round(cr,3)==round(nd.cr,3)
cross.diff[[i]] = cr[1,2]-nd.cr[1,2]
}
grads[[i]] = gradmi[1]-gradmi.nd[1]
}

#######
cd = unlist(cross.diff)
gd = unlist(grads)
hist(cd)
hist(gd)
hist(unlist(li.grads))
ncd = cd/min(cd)

us=unlist(sss)
plot(1:T,sss[[1]],ylim=(c(min(us),max(us))),type='l',col=0,xlab="time","ylab"="state")
for(i in 1:NC){
  thresh =1e-10
  lines(1:T,sss[[i]],col=(1+(abs(gd[i])>thresh)),lwd=(1+(abs(gd[i])>thresh)))
        #col=rgb(10*abs(ncd[i]),0,0))
}


# us=unlist(sss)
# plot(1:T,sss[[1]],ylim=(c(min(us),max(us))),type='l',col=0)
# for(i in 1:NC){
#   lines(1:T,sss[[i]],lwd=1)
#   #col=rgb(10*abs(ncd[i]),0,0))
# }
