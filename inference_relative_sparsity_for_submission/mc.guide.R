# runs coverage tables

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

# if use server set bh=1
bh=0 # need to set to 
print (c("bh=",bh))
if (bh){
args = commandArgs(trailingOnly=TRUE)
print(args)
# Rscript run.one.mc.R sel.tab.dir 8
res.dir = args[1]
param.path.sel = args[2]
}else{
  res.dir = exp.tag
  param.path.sel = param.path.sel
}

# check what selections were made in selection diagrams
vecs = readRDS(paste0(res.dir,"/unique.masks"))
# ignore all zero
vecs =vecs[unlist(lapply(vecs,function(x){!all(x==rep(0,length(vecs[[1]])))}))]
params.from.guide = readRDS(param.path.sel)

# get parameters given in run.mc.R
# some of these are live if not on server, so run through rstudio
print(c("params.from.guide",params.from.guide))
exp.tag = params.from.guide$exp.tag
param.path = paste0(exp.tag,"/param")

bh  = params.from.guide$bh

# run mc experiments for every selection
params = list()
for (ms in 1:length(vecs)){
  print(c("ms",ms))
  mysel = vecs[[ms]]
  params[[ms]] = list(deltas=params.from.guide$deltas,vecs=vecs,gammas=params.from.guide$gammas,
                      b0=params.from.guide$b0,n=as.integer(params.from.guide$n/2),tau=params.from.guide$tau,
                      T=params.from.guide$T,R=params.from.guide$R,
                      M=params.from.guide$M,K=params.from.guide$K,
                      sd.epsi=params.from.guide$sd.epsi,mean.first=params.from.guide$mean.first,
                      sd.first=params.from.guide$sd.first,
                      u=params.from.guide$u,plt.obj=FALSE,scale.s=params.from.guide$scale.s,
                      mysel=mysel,other.args=params.from.guide$other.args)
}

saveRDS(params,param.path)

# if on server, run in parallel
if (bh){
  dir.create('../../sbatch_scripts')
  dir.create(paste0(exp.tag,"/logs"))
  job.array.sbatch = "../../sbatch_scripts/run.para.sel.sbatch"
  custom.file = paste0(job.array.sbatch,'.',exp.tag)
  unlink(custom.file) 
  write(paste0("#!/bin/sh"),custom.file,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -n 1"),custom.file,append=TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs/log%a.txt"),custom.file,append=TRUE)
  write(paste0("#SBATCH -p standard --time=0-46:00:00"),custom.file,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file,append=TRUE)
  write(paste0("#SBATCH -a 1-",length(vecs)),custom.file,append=TRUE)
  write(paste0("module load r"),custom.file,append=TRUE)
  write(paste0("cd ",getwd()),custom.file,append=TRUE)
  write(paste0("Rscript ",getwd(),"/run.one.mc.R ",
               exp.tag," $SLURM_ARRAY_TASK_ID"),
        custom.file,append=TRUE)
  analysis.sbatch = "../../sbatch_scripts/analyze.mc.sbatch"
  
  custom.file2 = paste0(analysis.sbatch,'.',exp.tag)
  unlink(custom.file2) 
  write(paste0("#!/bin/sh"),custom.file2,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -p preempt --time=0-00:10:00"),custom.file2,append=TRUE)
  write(paste0("#SBATCH -n 1"),custom.file2,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file2,append = TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file2,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file2,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs/mc.table.txt"),custom.file2,append=TRUE)
  write(paste0("module load r"),custom.file2,append=TRUE)
  write(paste0("cd ",getwd()),custom.file2,append=TRUE)
  write(paste0("Rscript ",getwd(),"/analyze.mc.R ",
               exp.tag,"/sel.tables"),custom.file2,append=TRUE)
  
  pipeline =  paste0("../../sbatch_scripts/Infile.",exp.tag)
  unlink(pipeline)
  write(custom.file,pipeline,append=TRUE)
  write(custom.file2,pipeline,append=TRUE)
  system(paste("sbatch-pipeline",pipeline))
  
}else{
  sel.tables = list()
  for (ms in 1:length(vecs)){
    # all this will go in separate script
    param = readRDS(param.path)
    param.ms = param[[ms]]
    ftt2 = get.cov.table(param.ms)
    sel.tables[[ms]] = ftt2
    saveRDS(list(tab=ftt2,par=param[[ms]]),paste0(sel.tab.dir,"/",ms))
  }
  
  print(sel.tab.dir)
  fs = list.files(path=sel.tab.dir,full.names=TRUE)
  print(fs)
  fs = mixedsort(fs)
  exps = lapply(fs,readRDS) # should be in order
  
  for (el in exps){
    print(xtable(el$tab, 
                 caption = paste0("$n=$",el$par$n,", $T=$",el$par$T,", 
            $K=$",el$par$K, ", $M=$",el$par$M,
                                  ", Active=", paste0(el$par$mysel,sep=" ",collapse=""))),
          type="latex",sanitize.text.function = function(x){x},
          include.rownames=TRUE,digits=2)
  }
  
}
