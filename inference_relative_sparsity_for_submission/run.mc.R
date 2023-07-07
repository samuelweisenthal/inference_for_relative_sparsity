# this runs the simulation selection diagrams

library(latex2exp)
library(numDeriv)
library(Deriv)
library(xtable)
library(gtools)

source('get.var.R')
source('grad.R')
source('utils.R')
source('sim.R')
source('center.scale.R')
source('is.R') #vn
source('liops.R')
source('mc.utils.R')

# run on server
bh=0
#debug: set to 1 to run small experiment
debug =1
# Not scaling - explore more later
scale.s=0

#treatment effect
tau=.1

####
#other.args passed to functions; use.diff deprecated, z deprecated
# weighted var: use weighted IPW, dont.est.b = fixed behavior
# use.l is use KL penalty
# plot emp var if plot monte-carlo ci
# vmin = threshold. iff null, use num.se
# ci for ylim in coeff plots
# use.value.ci.for.ylim in value plots
# delta select ix
# gamma select ix
# C depreceated
other.args = list(use.diff=.1,z=0.2,weighted.var=1,dont.est.b=0,use.l=1,
                  plot.emp.var=1,vmin=NULL,num.se=2,use.ci.for.ylim=1,use.value.ci.for.ylim=1,
                  delta.select.hardcode=1,gamma.select.ix=1,C=1)
use.diff = other.args$use.diff
####

# sample size
n = 1000#400

#time steps
T=2#2

# MC iterations
M=5#00#00#20

# iterations for selection diagrams
plotsM = 3#3#0#0

# dimension of state
K=2#2

# lambda, gamma, delta grid
nlam =5
lambdas = seq(0.002,2,length.out=nlam)#seq(0.001,20,length.out=nlam)#seq(exp(log(0.001)),e$
#gammas = c(2,3,4)#c(1e-1)#seq(1,5,length.out=4)#c(1e-5,3.33e-2,6.67e-2,1e-1)#c(0.01,0.1,1)
gammas = seq(.01,2.6,length.out=3)[2]#c(1e-1)#seq(1,5,length.out=4)#c(1e-5,3.33e-2,6.67e-2,1e-1)#c(0.01,0.1,1)
deltas = c(0,1,2)[2]#c(0.5,1,2)

# which one to select: now passed above in other.args
lambda.select.ix = 2
gamma.select.ix = 2
gamma.select.ix = 2
lambda.select.ix = 2


if (debug){
# small experiment 
  n = 100#400
  T=2#2
  M=5#20
  plotsM = 5
  K=4#2#2
  nlam =4
  lambdas = seq(0.0001,1,length.out=nlam)#seq(0.001,20,length.out=nlam)#seq(exp(log(0.001)),e$
  other.args = list(use.diff=.1,z=0.2,weighted.var=1,dont.est.b=0,use.l=1,
                    plot.emp.var=1,vmin=NULL,num.se=.1,use.ci.for.ylim=1,use.value.ci.for.ylim=1,
                    delta.select.hardcode=1,gamma.select.ix=1,C=1)
  gammas = c(0,1,3)#c(1e-1)#seq(1,5,length.out=4)#c(1e-5,3.33e-2,6.67e-2,1e-1)#c(0.01,0.1,1)
  deltas = c(0.2,.3,3)
  lambda.select.ix = 1
  gamma.select.ix = 1
}

# put else sink
# Zou (adpaptive): pick a delta>0

# identifier of experiment
str.tag = "new"
exp.tag = paste0(str.tag,"Res.Dir.n=",n,"T=",T,"M=",M,"plotsM=",plotsM,"K=",K,
                 "gammaselix=",
                 gamma.select.ix,"lambda.select.ix=",lambda.select.ix,
                 "maxLam=",formatC(max(lambdas), format = "e", digits = 2),
                 "minLam=",formatC(min(lambdas),format="e",digits=2),
                 "maxGamma=",formatC(max(gammas), format = "e", digits = 2),
                 "minGamma=",formatC(min(gammas), format = "e", digits = 2),
                 "usediff=",use.diff,"tau=",tau)

# create directories
dir.create(exp.tag)
outfile = paste0(exp.tag,"/mainlog.txt")
plotfile = paste0(exp.tag,"/grid.mc.png")
resfile = paste0(exp.tag,"/outer.res")
param.path.sel = paste0(exp.tag,"/paramsel")
sel.tab.dir = paste0(exp.tag,"/sel.tables")
res.dir.sel = paste0(exp.tag,"/res.sel")
mylogdir = paste0(exp.tag,"/logs")
res.dir.sel.plots =  paste0(exp.tag,"/plots")
dir.create(res.dir.sel)
dir.create(res.dir.sel.plots)
dir.create(mylogdir)

if (!debug){
sink(file=outfile)
writeLines(readLines("run.mc.R"))
}

#sel.tables.path  = paste0(exp.tag,".sel.table")
dir.create(sel.tab.dir)

# define reward
#R = function(s,a,sp,T){
R = function(s,a,T){
  if (is.na(matrix(s)[1,])){
    r= 0
  }else{
    #r=-a*(s[2,] - 0.5*s[5,])
    r=-a*s[2,]   
    #r=-a*s
  }
  if(length(r)>1){
    print("r/state dimension mismatch.")
    browser()
  }
  r
}

u=1
b0=matrix(rep(0,K))
b0[1,]=-0.3
if (K>1){
b0[2,]=0.2
}

# define simulation settings
init.state.mean = matrix(rep(0,K))
init.state.var=diag(K)
trans.var = 1*diag(K)
sd.epsi=diag(rep(1,K))
mean.first=matrix(rep(0,K))
sd.first=diag(K)
#tau=2

rdigits =2
#eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T,tau=tau,R=R)

names = paste0("$S_{t,",1:K,"}$")

# put parameters into data structure
params.sel  =  list(deltas=deltas,gammas=gammas,lambdas=lambdas,
                    b0=b0,n=n,tau=tau,T=T,R=R,M=M,plotsM=plotsM,K=K,
                    sd.epsi=sd.epsi,mean.first=mean.first,
                    sd.first=sd.first,
                    u=u,plt.obj=FALSE,scale.s=scale.s,gamma.select.ix=gamma.select.ix,
                    lambda.select.ix=lambda.select.ix,exp.tag=exp.tag,
                    res.dir.sel = res.dir.sel,rdigits=rdigits,names=names,
                    #use.diff=use.diff,
                    bh=bh,
                    other.args=other.args)

saveRDS(params.sel,param.path.sel)

# if one server, call sbatch scripts to trun in parallel
if (bh==1){
  dir.create('../../sbatch_scripts')
  dir.create(paste0(exp.tag,"/logs.sel"))
  job.array.sbatch = "../../sbatch_scripts/run.sel.sbatch"
  custom.file.sel = paste0(job.array.sbatch,'.',exp.tag)
  unlink(custom.file.sel) 
  write(paste0("#!/bin/sh"),custom.file.sel,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -n 1"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -a 1-",plotsM),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs.sel/log%a.txt"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -p interactive --time=0-10:00:00"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.sel,append=TRUE)
  write(paste0("module load r"),custom.file.sel,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.sel,append=TRUE)
  write(paste0("Rscript ",getwd(),"/run.one.select.R ",
               param.path.sel," $SLURM_ARRAY_TASK_ID"),
        custom.file.sel,append=TRUE)
  #system(paste0("sbatch ",custom.file.sel))
  
  analysis.sbatch = "../../sbatch_scripts/analyze.sel.sbatch"
  
  custom.file.sel.2 = paste0(analysis.sbatch,'.',exp.tag)
  unlink(custom.file.sel.2) 
  write(paste0("#!/bin/sh"),custom.file.sel.2,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -p preempt --time=0-00:10:00"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH -n 1"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.sel.2,append = TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs/sel.table.txt"),custom.file.sel.2,append=TRUE)
  write(paste0("module load r"),custom.file.sel.2,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.sel.2,append=TRUE)
  write(paste0("Rscript ",getwd(),"/analyze.sel.R ", 
               exp.tag," ",
               res.dir.sel," ",res.dir.sel.plots),custom.file.sel.2,append=TRUE)
  
  mc.guide.sbatch = "../../sbatch_scripts/mc.guide.sbatch"
  
  
  custom.file.mc.guide = paste0(mc.guide.sbatch,'.',exp.tag)
  unlink(custom.file.mc.guide)
  write(paste0("#!/bin/sh"),custom.file.mc.guide,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -n 1"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/log.mc.guide.txt"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH -p interactive --time=0-10:00:00"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.mc.guide,append=TRUE)
  write(paste0("module load r"),custom.file.mc.guide,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.mc.guide,append=TRUE)
  write(paste0("Rscript ",getwd(),"/mc.guide.R ", 
               exp.tag," ", param.path.sel),custom.file.mc.guide,append=TRUE)
  
  pipeline =  paste0("../../sbatch_scripts/Infile.sel.",exp.tag)
  unlink(pipeline)
  write(custom.file.sel,pipeline,append=TRUE)
  write(custom.file.sel.2,pipeline,append=TRUE)
  write(custom.file.mc.guide,pipeline,append=TRUE)
  system(paste("sbatch-pipeline",pipeline))
  
}else{
  for (m in 1:plotsM){
    run.one.sel(params.sel,m)
  }
  analyze.sel.res(exp.tag,res.dir.sel,res.dir.sel.plots)
}
