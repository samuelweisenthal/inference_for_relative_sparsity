library(xtable)
source('mc.utils.R')
source('sim.R')
#sel.tables = list()
#for (ms in 1:length(vecs)){
  # all this will go in separate script
args = commandArgs(trailingOnly=TRUE)
print(args)
# Rscript run.one.mc.R sel.tab.dir 8
dire = args[1]
ms = as.integer(args[2])

param = readRDS(paste0(dire,"/param"))
print(param)
param.ms = param[[ms]]
ftt2 = get.cov.table(param.ms)
 # sel.tables[[ms]] = ftt2
saveRDS(list(tab=ftt2,par=param.ms),paste0(dire,"/sel.tables/",ms))
#}