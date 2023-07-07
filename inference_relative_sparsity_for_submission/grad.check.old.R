# check gradients against numerical
# note that actually m is not used anymore, but this does check all pieces of gradients that make up m,
# which are the pieces we need to check

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
#NC=1000
#for (i in 1:NC){
#  set.seed(i)
gamma=0.9
Ss=list(matrix(1),matrix(2))
As=list(1,0)
Rs=list(1,2)
T=length(Ss)
u=1
#### numerical dervis (better for multistage, not exact)
b=1
beta=2
sel=c(1)
sigmas=1
round(grad(kappa,x=beta,eta=b,Ss=Ss,As=As,sel=sel),7)==round(dthetakappa(beta,eta=b,Ss,As,sel),7)
(round(grad(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u,sigmas=sigmas,sel=sel),7)
  ==round(dbetami(beta,b,Ss,As,Rs,gamma,u=u,sigmas=sigmas,sel=sel),7))
(round(hessian(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,u=u,Rs=Rs,sigmas=sigmas,sel=sel),7)
  ==round(ddbetami(beta=beta,b=b,Ss=Ss,As=As,Rs=Rs,u=u,gamma=gamma,sigmas=sigmas,sel=sel),7))
round(dli(Ss,As,b,1),7)==round(grad(li,x=b,Ss=Ss,As=As,sel=1),7)
#round(ddlidbdb(1,1,1),7)==
round(ddli(Ss,As,b,1),3)==round(hessian(li,x=b,Ss=Ss,As=As,sel=1),3)
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
  if(K>1){
-s[2]*a
  }else{
    -s*a
  }
}

#eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T=50,tau=.1,R)
sel = c(0,1)#c(0,1)#c(0,1)#sample(c(0,1),K)
sel.b = rep(1,K)
sigmas=rep(1,K)
b0 = matrix(rep(0,K))
gamma=1
init.state.mean = matrix(rep(0,K))
init.state.var=.1*diag(K)
trans.var = .1*diag(K)
sd.epsi=diag(rep(1,K))
mean.first=matrix(rep(0,K))
sd.first=diag(rep(1,K))
tau=5#.1
T=2
#eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau=tau,R=R)

eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,
           sd.epsi=sd.epsi,mean.first=mean.first,sd.first=sd.first,seed=2)
# just generate

Ss=eps[[1]]$Ss
As=eps[[1]]$As
Rs=eps[[1]]$Rs

#sss[[i]] = unlist(Ss)#unlist.reshape(Ss)#[,2]
#T=100#length(Ss)
u=1
#### numerical dervis (better for multistage, not exact)
beta= matrix(rep(1,K))
sigmas=rep(1,K)
b=b0#matrix(c(-1))
round(grad(kappa,x=beta,eta=b,Ss=Ss,As=As,sel=sel),7)==round(dthetakappa(beta,b,Ss,As,sel),7)
gradmi.nd = grad(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u,sigmas=sigmas,sel=sel)
grad(vi,x=beta,b=b,Ss=Ss,As=As,Rs=Rs,u=u,sel=sel)
grad(li,x=beta,Ss=Ss,As=As,sel=sel)
gradmi = dbetami(beta,b,Ss,As,Rs,gamma,u,sigmas,sel)
print("dbetami")
round(gradmi.nd,7)==round(gradmi,7)

grad(kappa.b,x=b,Ss=Ss,As=As)
detakappa.b(eta=b,Ss=Ss,As=As)

print("dbri")
grad.b.ri.nd = grad(ri,x=b,beta=beta,Ss=Ss,As=As,sel=sel)
grad.b.ri=dbri(beta=beta,b=b,Ss=Ss,As=As,sel=sel)
round(grad.b.ri.nd,5) == round(grad.b.ri,5)

print("dbvi")
grad.vi.b.nd=grad(vi,x=b,beta=beta,Ss=Ss,As=As,Rs=Rs,u=u,sel=sel)
grad(vi.bfirst,x=b,beta=beta,Ss=Ss,As=As,Rs=Rs,u=u,sel=sel)
grad.vi.b=dbvi(beta=beta,b=b,Ss=Ss,As=As,Rs=Rs,u=u,sel=sel)
round(grad.vi.b.nd,5)==round(grad.vi.b,5)

print("ddbetami")
gradgradmi = ddbetami(beta=beta,b=b,Ss=Ss,As=As,Rs=Rs,gamma=gamma,u=u,sigmas=sigmas,sel=sel)
gradgradmi.nd = hessian(mi,x=beta,b=b,Ss=Ss,As=As,gamma=gamma,Rs=Rs,u=u,sigmas=sigmas,sel=sel)
(round(gradgradmi.nd,7)
  ==round(gradgradmi,7))
grad.li = dli(Ss,As,b,sel.b) 
grad.li.nd = grad(li,x=b,Ss=Ss,As=As,sel=sel.b)
#li.grads[[i]]=grad.li[1]-grad.li.nd[1]
round(grad.li,7)==round(grad.li.nd,7)
#round(ddlidbdb(1,1,1),7)==
gradgradli = ddli(Ss,As,b,sel.b)
round(gradgradli,4)==round(hessian(li,x=b,Ss=Ss,As=As,sel=sel.b),4)
#gnd.grad = function(b){grad(mi,x=beta,b=b,Ss=Ss,As=As,Rs=Rs,gamma=gamma)}
#nd.cross = grad(gnd.grad,x=beta) # deosn't work for vector inputs.
#round(nd.cross,5)==

# the cross products are more difficult to check,
# we will have to do it one at a time

# I think we assume
# d/db (d/dbeta_1 X d/dbeta_2 X)^T
# =((d/db_1 d/dbeta_1 X d/db_1 d/dbeta_2 X)^T, (d/db_2 d/dbeta_1 X d/db_2 d/dbeta_2 X)^T)
# so d/dbeta X is a column vector, and then d/db d/dbeta is a matrix with 
# b changing index over columns and beta over rows

if(K>1){
#cr=round(dbdbetami(beta,b,Ss,As,Rs,gamma,u,sigmas,sel),5)
cr=dbdbetami(beta,b,Ss,As,Rs,gamma,u,sigmas,sel)
  
beta1 = beta[1]
beta2 = beta[2]
b1 = b[1]
b2= b[2]

mi.ind.11 = function(beta1,b1,sigmas,sel){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u,sigmas,sel)
}
grad.m11 = function(b1,sigmas,sel){
  grad(mi.ind.11,x=beta1,b1=b1,sigmas=sigmas,sel=sel)
}

mi.ind.22 = function(beta2,b2,sigmas,sel){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u,sigmas,sel)
}

grad.m22 = function(b2,sigmas,sel){
  grad(mi.ind.22,x=beta2,b2=b2,sigmas=sigmas,sel=sel)
}

c22=grad(grad.m22,x=b2,sigmas=sigmas,sel=sel)

c11=grad(grad.m11,x=b1,sigmas=sigmas,sel=sel)

mi.ind.12 = function(beta1,b2,sigmas,sel){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u,sigmas,sel)
}

grad.m12 = function(b2,sigmas,sel){
  grad(mi.ind.12,x=beta1,b2=b2,sigmas=sigmas,sel=sel)
}

c12=grad(grad.m12,x=b2,sigmas=sigmas,sel=sel)

mi.ind.21 = function(beta2,b1,sigmas,sel){
  mi.ind(beta1,beta2,b1,b2,Ss,As,gamma,Rs,u,sigmas=sigmas,sel=sel)
}

grad.m21 = function(b1,sigmas,sel){
  grad(mi.ind.21,x=beta2,b1=b1,sigmas=sigmas,sel=sel)
}

c21=grad(grad.m21,x=b1,sigmas=sigmas,sel=sel)

# this seems to depend sometimes on the seed. if not equal then close
nd.cr=matrix(c(c11,c21,c12,c22),nrow=2,ncol=2)
print(sel)
print("Cross dbdbetami")
print(round(cr,3)==round(nd.cr,3))
print("it's not exact")
print(cr)
print(nd.cr)



}
