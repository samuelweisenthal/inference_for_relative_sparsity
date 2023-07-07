source('utils.R')
#source('mle.R')
source('is.R')
source('sim.R')
source('vops.R')
source('liops.R')
source('jointops.R')
source('expops.R')
source('cont.Ps.R')
source('plot.utils.R')
source('l1.R')
source('center.scale.R')

# runs one MC iteration 

args = commandArgs(trailingOnly=TRUE)
print(args)
save.direc = args[1]
my.mc.no = as.integer(args[2])

exp.list = readRDS(paste0(save.direc,"/param"))

told = Sys.time()
exp.res = wrap.jopt.lambda(exp.list[[my.mc.no]])

print("runtime")
rt=difftime(Sys.time(),told,units='mins')
print(rt)
save.direc
