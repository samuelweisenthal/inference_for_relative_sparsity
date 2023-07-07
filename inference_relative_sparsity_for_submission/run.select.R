library(xtable)
library(latex2exp)
source('mc.utils.R')
source('sim.R')

args = commandArgs(trailingOnly=TRUE)
print(args)

dire = args[1]

param = readRDS(dire)
exp.tag = param$exp.tag
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
plt.obj=param$plot.obj
scale.s=param$scale.s
mysel=param$mysel
gamma.select.ix = param$gamma.select.ix
lambda.select.ix = param$lambda.select.ix


eps=sim.db(b0=b0,N=n,tau=tau,T=T,R=R,
           sd.epsi=sd.epsi,mean.first=mean.first,sd.first=sd.first)

ps = check.pos(eps,paste0("S",1:K))

print(xtable(ps$dftb),
      type="latex",sanitize.text.function = function(x){x},
      include.rownames=TRUE,digits=2)

Ts=unlist(lapply(eps,function(x){length(x$As)}))
Ms=unlist.reshape(lapply(eps,function(x){x$Ms}))
#eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau=tau,R=R)
Ss = unlist.reshape.eps.ss(lapply(eps,function(x){x$Ss}))
rm = unlist.reshape(lapply(eps,function(x){x$Rs}))

nlam = length(lambdas)
ngam = length(gammas)
ndel = length(deltas)

plotfile = paste0(exp.tag,"/grid.mc.png")
resfile = paste0(exp.tag,"/outer.res")

selres = select.and.est(eps,b0,gammas,lambdas,deltas,names=paste0("S",1:K),
                     resfile=resfile,plotfile=plotfile,
                     scale.s=scale.s,lambda.select.ix=lambda.select.ix,
                     gamma.select.ix=gamma.select.ix,rdigits=2)



nplots = 2
res = 100*ngam
png(plotfile,height=500*ndel*nplots,width=800*ngam,res=res)
par(mfcol=c(ndel*nplots,ngam)) #vertical in gamma
#par(mfcol=c(ngam*5,ndel)) #horiz in gamma
# : bottom, left, top, and right.
par(mar=c(4.1,5.1,2,2))
ixs=plot.lambda.paths(selres$outer.res)
dev.off()

print("done")
saveRDS(selres$outer.res,resfile)


