library(latex2exp)
library(kableExtra)
library(xtable)
#library(neighbr)
library(predtools)
library(glmnet)
library(plotrix)
#library(reticulate) # didn't work. for reading python numpy arrays into r
source('utils.R')
source('mc.utils.R')
source('is.R')
source('mimic.utils.R')

# only micu (just a tag)
micu.only = 1

# if debug run small experiment
debug = 0
#don t change tag
tag = ""

# stages directory
stagesdirec = paste0("/stages",tag)
str.tag = "BehCon"

# which gamma and lambda to select if manual. Now given below in other.args
gamma.select.ix = 1
manu.lam.select = 3

####
# correspond to same args in run.mc.R
other.args = list(use.diff=0.01,z=0.2,weighted.var=1,dont.est.b=0,use.l=0,
                  plot.emp.var=0,vmin=NULL,num.se=.1,use.ci.for.ylim=0,use.value.ci.for.ylim=1,
                  delta.select.hardcode=1,gamma.select.ix=1,C=1)
use.diff = other.args$use.diff
####

max.lambda = 30#.02
nlam = 5 #10#10
lambdas= seq(1e-1,max.lambda,length.out=nlam)#exp(seq(log(2e-3),log(max.lambda),length.out=nlam))[2]#seq(0,max.lamda,length.out=nlam)
#lambdas = exp(seq(-50,log(max.lambda),length.out=nlam))
deltas = c(.5)#seq(1,2,length.out=3)#seq(.5,1.5,length.out=3)[3]
gammas = c(5)#exp(seq(-20,3,length.out=3))#exp(seq(-5,3,length.out=3))#exp(seq(-1.5,-0.5,length.out=3))#exp(seq(-1,0,length.out=3))#exp(seq(-22,19,length.out=3))#exp(seq(-2,-1,length.out=3))#exp(seq(-2,0,length.out=3))#exp(seq(-2,0,length.out=3))[1]#exp(seq(-2,1,length.out=3))#exp(seq(-10,0,length.out=3))#exp(seq(-26,2,length.out=3))

if (debug){
  nlam = 2#5#10 #10#10
  lambdas= seq(0.1,60,length.out=nlam)#exp(seq(-10,log(max.lambda),length.out=nlam))#seq(0,max.lamda,length.out=nlam)
  #lambdas = exp(seq(-50,log(max.lambda),length.out=nlam))
  deltas = seq(.5,1.5,length.out=3)
  gammas = exp(seq(-26,2,length.out=3))
  
}

#gammas = c(1e-50)
nlam = length(lambdas)
ngam = length(gammas)
ndel = length(deltas)

# how many covariates
ncov = 9#9#9 # 8 is max
# what we want to focus on. Makes a two stage decision. 3 is actually beginning (there are buffer stages that are like a few seconds)
start.t=3 # beginning of traj
end.stage=4 # end of traj. we end early. so this is like first 30 min after hyp

# tag
exp.tag = paste0("tag=",str.tag,"usediff=",use.diff,",start.t",start.t,"endStage",
                 end.stage,"gammaselix=",
                 gamma.select.ix,"manu.sel=",manu.lam.select,"nlam=",
                 nlam,"minlam=",min(lambdas),"maxlam=",max(lambdas),"gammas=",
                 paste(formatC(gammas,digits=2,format="e"), collapse = '.'),
                 "ncov=",ncov,"use.dff=",use.diff)

careunit = 'MICU'
print("Need to make sure in proper directory!!!")
if(getwd()=="/Volumes/projects/Latent/sam/thesis_work/relative_sparsity/code_betareg/mdp_betareg_KL"){
  
  if (micu.only){
  pth='/Users/ANONYMIZED/Box/MIMIC/mimic-iii-v1-4/hypotension-RL/model-data2'  
  #pth='/Users/ANONYMIZED/Box/MIMIC/mimic-iii-v1-4/hypotension-RL/model-data3'
  }else{
  
  pth='/Users/ANONYMIZED/Box/MIMIC/mimic-iii-v1-4/hypotension-RL/model-data-allicu'
  }
}else{
  print("stop")
  browser()
  pth='/scratch/ANONYMIZED/mimic2'
  
  sink(file = paste0(exp.tag,".mimic.txt"), append = FALSE, 
       type = c("output", "message"),split = FALSE)
}



writeLines(readLines("mimic.R"))

# unit test, to make sure correct info for particular patient
 id.check = 52
# s[id.check,]
# a[id.check]
# r[id.check]

# get stage info
filenames <- list.files(paste0(pth,stagesdirec), pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)

png("first.maps.png")
hist(ldf[[1]][,'map'],xlab="first maps")
dev.off()

# get covariates
ldf = lapply(ldf,as.matrix)
# what is available (change in np_load.py)
nStages=length(ldf)
#first.map.zero=(
#  (ldf[[1]][,'first_lowMAPix']==0) 
#  & (!is.na(ldf[[1]][,'first_lowMAPix'])))
state.dim = dim(ldf[[1]])[2] #K
nObs = dim(ldf[[1]])[1] #N
#data.array = array(NA,c(nObs,nStages,state.dim)) #T

# get rewards, actions, and times
msr =  read.csv(paste0(pth,'/',tag,'multistage_rewards.csv'),header=TRUE, row.names = 1)
msa =  read.csv(paste0(pth,'/',tag,'multistage_actions.csv'),header=TRUE, row.names = 1)
mst =  read.csv(paste0(pth,'/',tag,'multistage_times.csv'),header=TRUE, row.names = 1)
#msva = read.csv(paste0(pth,'/icu_startmultistage_vasoamts.csv'),header=TRUE, row.names = 1)

# 
# unit test
msr[id.check,]
msa[id.check,]
#msva[id.check,]
mst[id.check,] # note that these times correspond to the rewards and actions, not states
c(ldf[[1]][id.check,]['map'],
  ldf[[2]][id.check,]['map'],ldf[[3]][id.check,]['map'])
### vaso effect from time 1 to 2
init.map = ldf[[1]][,'map']
final.map = ldf[[2]][,'map']
vaso=msa[,1]
# some exploratory diagrams
df = as.data.frame(cbind(ldf[[1]],vaso,final.map))
#summary(lm(final.map~.,data=df))



no.first = is.na(msa[,1])
print("% without first (why are these occuring?)")
sum(no.first)/length(no.first)*100


misa = rep(0,length(msa[,1]))

for (i in start.t:(end.stage)){
  print("stage")
  print(c("n missing",sum(is.na(msa[,i]))))
  misa = misa + is.na(msa[,i])
}


#ms=1
#if(ms){
is.misa=misa>0
print("subsetting by first map zero")
take = !is.misa #& first.map.zero


print(c("total missing",sum(1-take)))
print(c("n",sum(take)))
ixs = 1:dim(msa)[1]
newixs = ixs[take]


for (i in 1:length(ldf)){
  print(dim(ldf[[i]]))
}

# just take decisions betwwen start.t and end.stage
# starting at time start.t
print("subsetting by first map zero")
msr = msr[,start.t:(end.stage),drop=FALSE]
msa = msa[,start.t:(end.stage),drop=FALSE]
ldf = ldf[start.t:(end.stage)]
mst = mst[,start.t:(end.stage),drop=FALSE]

msr=msr[take,,drop=FALSE]
msa = msa[take,,drop=FALSE]
mst = mst[take,,drop=FALSE]
ldf=lapply(ldf,function(x){x[take,,drop=FALSE]})

png("stagesActRew.png",width=1000,height=dim(msr)[2]*1000,res=150)
par(mfrow=c(dim(msr)[2],1))
for (i in 1:dim(msr)[2]){
  boxplot(msr[,i]~msa[,i],main=paste("stage",i),xlab="Action",ylab="Reward")
}
dev.off()


#}

#unit test again
id.check = which(newixs==id.check)

msr[id.check,]
msa[id.check,]
#msva[id.check,]
mst[id.check,] # note that these times correspond to the rewards and actions, not states


png("act.rew.over.time.png",width=1000,height=1000,res=150)
avact = apply(msa,2,mean) 
avre = apply(msr,2,mean)
par(mfrow=c(2,1))
plot(c(mst[1,]),avact,xlab="time (hrs)",ylab="prop. given vasos")
plot(c(mst[1,]),avre,xlab="time (hrs)",ylab="avg reward")
dev.off()




# Excluding categorical variables for now
#I also noticed that Joe’s paper only uses 9 variables, so maybe we could follow him. 
#He uses MAP, heart rate, urine output, lactate, Glasgow coma score, 
# serum creatinine, FiO2, total bilirubin, and platelets count. 







# Futoma et al. 2020 variables
cov.of.int = c("map","hr","urine","lactate","GCS","creatinine","fio2",
               "bilirubin_total","platelets")[1:ncov]


ldf2=lapply(ldf,function(x){x[,cov.of.int,drop=FALSE]})

cov.of.int  =  c("MAP","HR","urine","lactate","GCS","creatinine","Fio2",
                 "bilirubin","platelets")[1:ncov]
#colnames(s) = cov.of.int
for (i in 1:length(ldf2)){
  colnames(ldf2[[i]])= cov.of.int
  
}

# turn raw data into episodes
mseps = get.eps.ms(ldf2,msa,msr)

# check for positivity violations
check.pos(mseps,cov.of.int)

mc=check.pos.cal(mseps,cov.of.int,pen.b=0,train.test=1,plot.=TRUE,usesampler=FALSE,other.args)

if (micu.only==0){
calname = "calc.png"
}else{
  calname="calcMICU.png"
}

# make calibration curve on resampled data and plot
plot.curve=1
if (plot.curve){
  png("calcMICU.png",res=150,height=800,width=800)
  
  M=100
  txs = matrix(NA,nrow=M,ncol=2)
  test.cal.plots = list()
  for (i in 1:M){
    mr = check.pos.cal(mseps,cov.of.int,pen.b=0,plot.=0,usesampler=TRUE,other.args=other.args)
    txs[i,]=c(mr$IPTW.treat.eff.train,mr$IPTW.treat.eff.test)
    test.cal.plots[[i]] = mr$p.test
  }
  md=matrix(apply(txs,2,mean),nrow=1,ncol=2)
  colnames(md) = c("IPTW treat. eff. (train)","IPTW treat. eff.(test)")
  kable(md)
  
  av.cal.test=Reduce('+',test.cal.plots)/M
  #plot(av.cal.test$predRate,av.cal.test$obsRate,xlim=c(0,0.4),ylim=c(0,0.4))
  #lines(av.cal.test$predRate,av.cal.test$obsRate_UCL)
  #lines(av.cal.test$predRate,av.cal.test$obsRate_LCL)
  #
  plot(c(0,1),c(0,1),type='l',xlim=c(0,0.4),ylim=c(0,0.4),
       xlab=TeX("Estimated ($\\pi_{b_n}(A_0=1|S_0=s_0)$)"),
       ylab=TeX("Observed"), main=paste0("Average Calibration Curve Over ",M, " Resamples"))
  plotCI(av.cal.test$predRate, av.cal.test$obsRate, ui=av.cal.test$obsRate_UCL, li=av.cal.test$obsRate_LCL,add=TRUE)
  dev.off()
}


# if debug just look at subset
if (debug){
  eps.mimic = mseps[1:1000]
}else{
  eps.mimic=mseps #mseps#[1:1000]#halfds
}

# scale covariates (don't center, old name)
sc=center.scale.s(eps.mimic)
smm=sc$sm.cs

# look at scaled variables
png("ScaledCovHist.png",width=1000,height=length(cov.of.int)*1000,res=150)
par(mfrow=c(length(cov.of.int),1))
for (i in 1:length(cov.of.int)){
  hist(smm[,i],main=paste(cov.of.int[i]))
}
dev.off()

print("SET Sigmas to 1")

ngam*ndel*nlam
scale.s = 1 # for real data, must scale
# for mimic, b0 doesn't matter. only for coverage results

# run analysis
ixs=select.and.est(eps.mimic,b0=rep(1,dim(eps.mimic[[1]]$Ss[[1]])[1]),
                   gammas=gammas,lambdas=lambdas,
                   deltas=deltas,names=cov.of.int,
                   resfile=paste0(exp.tag,"outer.res.mimic"),
                   plotfile=paste0(exp.tag,"grid.mimic.png"),
                   scale.s=scale.s,
                   lambda.select.ix=manu.lam.select,
                   gamma.select.ix=gamma.select.ix,
                   other.args=other.args)
print(ixs)
