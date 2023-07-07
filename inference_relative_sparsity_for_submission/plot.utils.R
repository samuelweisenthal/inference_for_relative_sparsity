# this is the workhorse of analysis.R
# it contains functions for plotting and some analysis
# it is somewhat repetitive
# it is organized by plot type

plot.17 = 0 # whether to highlight 17, necessary for Appendix experiment 2

get.summ.stat = function(ys,getpm,xs){
  mean.mean.table.s.bar = tapply(ys,xs,getpm)#mean)
  my.lgrid = as.numeric(names(mean.mean.table.s.bar))
  list(summ.stat=mean.mean.table.s.bar,lgrid=my.lgrid)
}

wilcox.test.wrapper = function(v,conf.int,conf.level,exact){

  op=wilcox.test(v,conf.int,conf.level,exact)
  
  list(conf.int=op$cont.int,estimate=op$estimate)
}

get.q.ci = function(ys,xs,cl=0.95){
  #browser()
  wc = tapply(ys,xs,wilcox.test,conf.int=TRUE,conf.level=cl,exact=FALSE)
  lwc =unlist(lapply(wc,function(x){x$conf.int[1]}))
  uwc =unlist(lapply(wc,function(x){x$conf.int[2]}))
  upper.95 = tapply(ys,xs,quantile,probs=c(.975))
  lower.95 = tapply(ys,xs,quantile,probs=c(.025))
  estimate = unlist(lapply(wc,function(x){x$estimate}))
  list(lwc=lwc,uwc=uwc,upper.95=upper.95,lower.95=lower.95,estimate=estimate)
}

getc = function(x,cl){
  wilcox.test(x,conf.int = TRUE,conf.level = cl,exact=FALSE)$conf.int
  }
getpm = function(x){
  wilcox.test(x,conf.int = TRUE,exact=FALSE)$estimate
  }


scale = function(x,xs,ys){
  slope = ((min(xs)-max(xs))/(min(ys)-max(ys)))
  s.y = (x-(min(xs)-min(ys)*slope))/slope
  s.y
}

lamb.plot=function(theta.stars,parbs,lambdas,HAL=FALSE){

  # bottom left top right
  par(mar=c(4,6,1,1))
  cex.lab.size = 1.7
  par(mfrow=c(1,1))
  theta.stars.m = unlist.reshape(theta.stars)
  parbs.m = unlist.reshape(parbs)
  sds = unlist(lapply(exps,function(x){x$input.ep.seed}))
  ####
  df = cbind(theta.stars.m-parbs.m)

  df[sds==1,]
  lm.u.s = sort(unique(lambdas))
  get.first = function(x,lm.u.s){

    for (i in 1:length(x)){
      first=i
      if(is.na(abs(x[i]))){
        print("xi is na")
        browser()
      }
      if (abs(x[i])<1e-3){
        break
      }
      
    }
    lm.u.s[first]
  }
  firsts=list()
  for (j in 1:length(unique(sds))){
    firsts[[j]] = apply(df[sds==j,],2,get.first,lm.u.s=lm.u.s)
  }

  fm = matrix(unlist(firsts),nrow=length(unique(sds)),ncol=dim.s$optim.exp,byrow=TRUE)
  #hist(fm[,1],xlim=c(min(fm),max(fm)),col=1,xlab=TeX("$\\lambda_{i}^*$"),main="")
  cex.lab.size = 1.7

  ylimt=3.5#max(density(fm[,2])$y,from=min(lambdas),to=max(lambdas))

  densities = list()
  for (i in 1:dim.s$optim.exp){
    densities[[i]] = density(fm[,i],from=min(lambdas),to=max(lambdas))$y
  }
  myde = unlist(densities)
  ylimt = max(log(myde+1e-20))
  fmd = density(fm[,2],from=min(lambdas),to=max(lambdas))
  
  xlim=c(fmd$x[2],max(lambdas))
  xlim = c(min(lambdas),max(lambdas))
  fmd = density(fm[,1],from=min(lambdas),to=max(lambdas))

  #if(HAL){
    logv = 'x'#'xy'
    lty=1
    ylab="Log density"
    ymin  = min(log(myde+1e-20))
  #}else{
  #  logv='x'
  #  lty=2
  #  ylab="Density"
  #  ymin = 0
  #}
  #ymin=-1
  alpha=0.05
  fmd$y = log(fmd$y+1e-20)
  plot(fmd$x,(fmd$y),ylim=c(ymin,ylimt),xlim=xlim,
       col=rgb(red=0,green=0,blue = 0, alpha=alpha),ylab=ylab,
       xlab=TeX("$\\zeta_{k}^*$"),main="",type='l',cex.lab=cex.lab.size,
       cex=1,lty=lty,log=logv)
  
  if(HAL){
    
    for (i in 1:(sort(dim.s$optim.exp))){

      #hist(fm[,i],add=TRUE,col=i,xlab="")
      fmd = density(fm[,i],from=min(lambdas),to=max(lambdas))
      
      if (sum(fmd$y<0)>0){
        browser()
      } 
      fmd$y = log(fmd$y+1e-20)
      if (i>(dim.s$optim.exp/2)){
        col = rgb(red = 1, green = 0, blue = 0, alpha=alpha)
      }else{
        col = rgb(red=0,green=0,blue = 0, alpha=alpha)
      }

      lines(fmd$x,fmd$y,lty=1,col=col)
    }
    legend("top",legend=c(expression(paste( k, " > ", N)),
                               expression(paste( k <= N))
    ),
    lty=c(1,1),
    col=c(2,1),lwd=c(1,1))

  }else{
  if (dim.s$optim.exp>2){
  for (i in 3:(dim.s$optim.exp)){
    print(i)
    #hist(fm[,i],add=TRUE,col=i,xlab="")
    fmd = density(fm[,i],from=min(lambdas),to=max(lambdas))
    fmd$y = log(fmd$y+1e-20)
    lines(fmd$x,fmd$y,lty=2)
  }
  }
  fmd = density(fm[,2],from=min(lambdas),to=max(lambdas))
  fmd$y = log(fmd$y+1e-20)
  lines(fmd$x,fmd$y,lty=1,col=2,lwd=2)
  #legend("topleft",legend=paste0("i=",1:dim.s$optim.exp),lty=1:dim.s$optim.exp,
  #       col=1:dim.s$optim.exp)
  legend("bottom",legend=c("k = 2",expression("k"!= 2)),lty=c(1,2),
         col=c(2,1),lwd=c(2,1))
  }
  #legend("topright",legend=c(expression(SSE[V[lambda]]),
  #                           expression(paste("|","!=","*",-V[lambda],"|"))),lty=1:2)
  
}



lasso.plot=function(theta.stars,parbs,lambdas,opfile,HAL=FALSE,BAS=FALSE){

  # bottom left top right
  par(mar=c(4,6,1,1))
  cex.lab.size = 1.7
  par(mfrow=c(1,1))

  theta.stars.m = unlist.reshape(theta.stars)
  parbs.m = unlist.reshape(parbs)
  
  df = (
    #abs(
    theta.stars.m-parbs.m
    #)+1e-100
    )

  meds = apply(df,2,get.summ.stat,getpm,lambdas)
  dim.s = length(theta.stars[[1]])
  if(BAS){
    dim.s.raw=dim.s/BAS
    cols=rep(1:dim.s.raw,each=BAS)
  }
  par.in = list()
  for (i in 1:dim.s){
    par.in[[i]] = get.q.ci(df[,i],lambdas)
  }
  ####
  
  ####

  ylim = c(min(df),max(df))
  # the actual point estimates have so much variance, use quantiles to get bounds of plot?
  # I need to discuss about this. we don't want to "hide" outliers, but
  # there is also clearly a pattern here that we want to show..
  lwcs = unlist(lapply(par.in,function(x){x$lwc}))
  uwcs = unlist(lapply(par.in,function(x){x$uwc}))
  
  # remove NA bc sometimes there is no variation in numbers, 
  # need finite bounds for ylim of plot
  # (eg if diff btwn beta and b always 0)
  # so cannot construct CI
  # this depends on the dimension of the state, strength of penalty
  # (eg lgrid)
  # and values of the behavioral policy
    
  ylim2 = c(min(lwcs,na.rm=TRUE),
  max(uwcs,na.rm=TRUE))
#  if (is.na(ylim2[1])){
#    print("na.ylim")
# }
  
  diff.med = apply(df,2,get.summ.stat,getpm,lambdas)
  lgrid = diff.med[[1]]$lgrid
  xlim=c(min(lambdas),max(lambdas))
  #ok to plot 0 0 like this?
  ylab=TeX("$(\\beta_{n,\\gamma,\\lambda})_k-(b_{n})_{k}$")
  if(HAL){
    d=unlist(lapply(diff.med,function(x){x$summ.stat}))
    ylim2=c(min(d,na.rm=TRUE),max(d,na.rm=TRUE))
    ylab=TeX("$\\hat{\\phi}_{k,h}-\\hat{\\eta}_{k,h}$")
  }else if (BAS){
    ylab=TeX("$\\hat{\\phi}_{k,h}-\\hat{\\eta}_{k,h}$")
  }
  plot(lambdas,df[,2],ylim=ylim2,xlim=xlim,col="white",
       xlab = TeX("$\\lambda$"),ylab=ylab,cex.lab=cex.lab.size,
       #log='xy',
       log='x',
       cex=1,panel.first = c(abline(h = 0,lty=1)))
  
  if(plot.17){
    legend("topright",legend=c(expression("k = 2"),
              
                               expression(paste("k" != 2,", 17")),
                               expression("k = 17")
    ),
    lty=c(1,2,4),
    col=c(2,1,17),lwd=c(3,1,3))
  }else{
   if(HAL){
     legend("topright",legend=c(expression(paste("", h, " > ", N)),
                                expression(paste("", h <= N))
     ),
     lty=c(1,1),
     col=c(2,1),lwd=c(1,1))
     
   }else if(BAS){
     legend("topright",legend=paste("k=",1:dim.s.raw),
     lty=rep(1,dim.s.raw),
     col=1:dim.s.raw,lwd=rep(1,dim.s.raw))
   }
   else{
  legend("left",legend=c(expression("Pmedian k = 2"),
                                                       expression("95% CI, k = 2"),
                                                           expression("Pmedian k" != 2)
                                ),
           lty=c(1,1,2),
            col=c(2,2,1),lwd=c(3,1,0.5))

  
   }
  }
  

  if(plot.17){
  lines(lgrid,diff.med[[17]]$summ.stat,lty=4,cex=1,col=17,lwd=3)
  #lines(lgrid,par.in[[17]]$uwc,lty=2,cex=1,col=17,lwd=1)
  #lines(lgrid,par.in[[17]]$lwc,lty=2,cex=1,col=17,lwd=1)
  }
  if (HAL){
    se=sort(1:dim.s,decreasing=TRUE)
  }else{
  se = 1:dim.s}
  for(j in se){
      #if (j!=2){
      if (plot.17){
      if (j==17){}else{
        lines(lgrid,diff.med[[j]]$summ.stat,lty=2,cex=1,lwd=.75)#,log='x')
      }
      }else{
      
        if(HAL){
 
          col=((j>(dim.s/2))+1)
          lty=1
        }else if(BAS){
          col = cols[j]
   
          lty=1
        }else{
          col=1
          lty=2
        }
        
      lines(lgrid,diff.med[[j]]$summ.stat,lty=lty,cex=1,lwd=.75,col=col)#,log='x')
      }
      #}
  }
  if(!HAL&&!BAS){
  lines(lgrid,diff.med[[2]]$summ.stat,lty=1,cex=1,col=2,lwd=3)
  if (!plot.17){
  lines(lgrid,par.in[[2]]$uwc,lty=1,cex=1,col=2,lwd=1)
  lines(lgrid,par.in[[2]]$lwc,lty=1,cex=1,col=2,lwd=1)
  }
  }
  ####################3
  #"95CI for B2: ",
  nfCI = 0
  if (nfCI){
  nd=2
  par.in = list()
  for (i in 1:dim.s){
    par.in[[i]] = get.q.ci(df[,i],lambdas,cl=(1-0.05/dim.s))
  }
  
  df = cbind(round(par.in[[2]]$estimate,nd),
        paste0("(",round(par.in[[2]]$lwc,nd),
               ",",round(par.in[[2]]$uwc,nd),")"))
  
  myci = paste0(round(par.in[[2]]$estimate,nd),
                "(",round(par.in[[2]]$lwc,nd),",",round(par.in[[2]]$uwc,nd),")")
  myci = cbind(round(par.in[[2]]$estimate,nd),round(par.in[[2]]$lwc,nd),round(par.in[[2]]$uwc,nd))
  write.csv(paste0(save.direc,"/CIs2.txt"),myci)
  }
}

plot_=function(theta.star.value.train,theta.star.like.train,
               train.test=FALSE,beta.plots=FALSE,title=FALSE){
  #beta.plots=FALSE
  #if(beta.plots){fi=4}else{fi=2}
  # bottom left top right
  
  #par(mfrow=c(fi,2))
  par(mfrow=c(2,1))
  cex.lab.size = 1.7
  if (disc.a){
    doty = TeX("$\\hat{\\pi}(A=1|S=\\bar{s})$")
    doty.beta1 =  TeX("$\\beta_{n,\\lambda,1}$")
    doty.beta2 =  TeX("$\\beta_{n,\\lambda,2}$")
    doty.beta.val =  TeX("$\\hat{V}_{\\pi_{\\hat{\\theta}}}$")
    vy = TeX("$SSE_{\\hat{\\pi}(A=1|S=\\bar{s})}$") 
    
  }else{
    doty = TeX("$\\pi_{\\hat{\\theta}}(A=\\bar{A}|S=\\bar{S})$")
    doty.beta =  TeX("$\\beta_{n,\\gamma,\\lambda}$")
    doty.sigma =  TeX("$\\hat{\\sigma}$")
    vy = TeX("$SSE_{\\pi_{\\hat{\\theta}}(A=\\bar{A}|S=\\bar{S})}$")    
    vy.sigma = TeX("$SSE_{\\hat{\\sigma}}$")
    
  }
  vy.beta = TeX("$SSE_{\\beta_{n,\\gamma,\\lambda}}$")

  
  doty.theta.val.tr =  TeX("$V_n(\\beta_{n,\\gamma,\\lambda},b_n)$")
  doty.theta.like.tr =  TeX("$\\tilde{L}_n(\\beta_{n,\\gamma,\\lambda})$")

  
  
  vy.theta.val.tr = TeX("$SSE_{V_n(\\beta_{n,\\gamma,\\lambda},b_n)}$")
  vy.theta.like.tr = TeX("$SSE_{\\tilde{L}_n(\\beta_{n,\\gamma,\\lambda})}$")

  
  
  xlab = TeX("$\\lambda$")
  
  # need legend eventually.
  plot.dots = function(xs,ys,xlab,ylab,leg){
    # bottom left top right
    if (!leg){par(mar=c(4,6,1,1))}else{par(mar=c(1,6,4,1))}
    ob = get.q.ci(ys,xs)
    upper.95=ob$upper.95
    lower.95=ob$lower.95
    lwc = ob$lwc
    uwc = ob$uwc
    if(is.na(lwc[1])){
      print("Lower Wilson CI is NA. Maybe because all tied or all zero")
      print("Won't plot")
      browser()
    }
    alpha=.1
    if (dim.s$raw==1){conv=0}
    col.v = c(rgb(red=0,green=0,blue = 0, alpha=alpha),rgb(red=1,green=0,blue = 0, alpha=alpha))[conv+1]
    #plot(1e-50,1e-50, col='white',#
    loc = "bottomleft"  
    if(leg){
      title=title
      if(title=='Test'){loc = "bottomleft"}
      if(title=='Train'){loc= "bottom"}
      
      }else{title=NULL}
    plot(xs,ys,col=col.v,
          xlab=xlab,
         ylab=ylab, xlim=c(min(lambdas),max(lambdas)),ylim=c(min(lwc),max(uwc)),
         # I wish I could plot this on log scale y. should have chosen problem with nonneg value.
         main=title,cex=1,cex.lab=cex.lab.size,log='x') # too much info

    
    if(leg){legend(loc,legend=c("Pmedian","95% CI",expression(lambda*"*")),lty=c(1,2,4),
                   col=c(1,2,1),lwd=c(3,2,3))}
    nvrep = length(lambdas)/length(unique(lambdas))
    mean.mean.table.s.bar = tapply(ys,xs,getpm)#mean)
    #upper.95 = tapply(ys,xs,quantile,probs=c(.975))
    #lower.95 = tapply(ys,xs,quantile,probs=c(.025))
    sds = tapply(ys,xs,sd)

    real.mean = tapply(ys,xs,mean)#mean)
    zt = qt(0.975,df=(nvrep-1))
    # manual CI. then decided to use wilcox.text
    low  = real.mean + zt*(sds/sqrt(nvrep))
    high = real.mean - zt*(sds/sqrt(nvrep))
    
    
    my.lgrid = as.numeric(names(mean.mean.table.s.bar))
    mi=NULL
    if (leg){
      tomax = mean.mean.table.s.bar - (mean.mean.table.s.bar-lwc)
      mi = which(tomax==max(tomax))
      
      #mtext(expression(lambda[n]),side=2,at=my.lgrid[mi])
      abline(v=my.lgrid[mi],lty=4,col=1,lwd=3)
      #abline(v=ls[ix],lty=3)
      
    }
    
    #lines(my.lgrid,low,col=3)
    #lines(my.lgrid,high,col=3)
    lines(my.lgrid,lwc,lty=2,col=2,lwd=2)
    lines(my.lgrid,uwc,lty=2,col=2,lwd=2)
    lines(my.lgrid,mean.mean.table.s.bar,col=1,lty=1,lwd=3)#,log='x')
    
    #loc = "bottomleft"
    
   mi 
  }
  
  plot.dots.train.test = function(xs,ys.tr,ys.te,xlab,ylab){
    col.v = c("grey","red")[conv+1]
    plot(xs,ys.tr, xlab=xlab,
         ylab=ylab,
         main="",cex=1,cex.lab=cex.lab.size,col=col.v,log='x') # too much info
    points(xs,ys.te,col='cadetblue1')#,log='x')
    
    tr = get.summ.stat(ys.tr,getpm,xs)
    lines(tr$lgrid,tr$summ.stat,lty=1)#,log='x')
    te = get.summ.stat(ys.te,getpm,xs)
    lines(te$lgrid,te$summ.stat,lty=2)#,log='x')
    
    
  }
  
  plot.vars = function(xs,ys,xlab,ylab,leg){
    
    o = get.summ.stat(ys,function(x){sqrt(var(x))},xs)
    mean.table.s.bar = o$summ.stat#tapply(ys,xs,function(x){sqrt(var(x))})
    #my.lgrid = as.numeric(names(mean.table.s.min))
    my.lgrid = o$lgrid#as.numeric(names(mean.table.s.bar))
    if(leg){

    alltog = c(mean.table.s.bar)
    ylim = c(min(alltog),max(alltog))


    }else{
      ylim=NULL
    }

    plot.vars=FALSE
    if (plot.vars){
    plot(my.lgrid,mean.table.s.bar, xlab=xlab,
         ylab=ylab,ylim=ylim,
         main="",cex=1,cex.lab=cex.lab.size,type='l',log='xy') 
    if (leg){
    legend("bottomleft",legend=c("SSE (MC)","E[(V-barV)^2] (1 Dataset)"),pch=c(NA,1),
           lty=c(1,NA))
    points(lambdas,sd.vs)
    }
    }
    
  }
  
  
  plot.dots.and.vars = function(ys,dot.lab,var.lab,leg=FALSE){
    mi=plot.dots(xs=lambdas,ys=ys,xlab=xlab,ylab=dot.lab,leg)
    plot.vars(xs=lambdas,ys=ys,xlab=xlab,ylab=var.lab,leg)
    mi
  }
  
  
  if(train.test){

    par(mfrow=c(1,2))
    plot.dots.train.test(xs=lambdas,ys.tr=theta.star.like.train,
                         ys.te=theta.star.like.test,xlab="",ylab="")
    plot.dots.train.test(xs=lambdas,ys.tr=theta.star.value.train,
                         ys.te=theta.star.value.test,xlab="",ylab="")
  }else{

  mi = plot.dots.and.vars(ys=theta.star.value.train,
                     dot.lab=doty.theta.val.tr,var.lab=vy.theta.val.tr,leg=TRUE)
  
  #plot.dots.and.vars(ys=theta.star.value.test,
  #                   dot.lab=doty.theta.val.te,var.lab=vy.theta.val.te)
  
  
  plot.dots.and.vars(ys=theta.star.like.train,
                     dot.lab=doty.theta.like.tr,var.lab=vy.theta.like.tr,leg=FALSE)

  #if(disc.a){
    #print("I am plotting hardcoded 2 dim betas! fix in plot.utils.R")
  # if(beta.plots){  
  # plot.dots.and.vars(ys=beta.stars1,dot.lab=doty.beta1,var.lab=vy.beta,leg=FALSE)
  #   plot.dots.and.vars(ys=beta.stars2,dot.lab=doty.beta2,var.lab=vy.beta,leg=FALSE)
  # }
  #}else{
    #plot.dots.and.vars(ys=beta.stars,dot.lab=doty.beta,var.lab=vy.beta,leg=FALSE)
  #}
  #plot.dots.and.vars(ys=theta.star.like.test,
  #                   dot.lab=doty.theta.like.te,var.lab=vy.theta.like.te)
  }
  
  if (beta.plots){
    par(mfrow=c(2,1))
    if(disc.a){
      print("I am plotting hardcoded 1 dim betas! fix in plot.utils.R")
      plot.dots.and.vars(ys=beta.stars1,dot.lab=doty.beta1,var.lab=vy.beta)
    
      #plot.dots.and.vars(ys=beta.stars2,dot.lab=doty.beta2,var.lab=vy.beta)
    }else{
      plot.dots.and.vars(ys=beta.stars,dot.lab=doty.beta,var.lab=vy.beta)
    }
  }
  mi
}


lasso.plot2=function(theta.stars,lambdas,opfile){
  # plots lasso plot with absolute values. sort of nice
  par(mar=c(4,3,1,1))
  cex.lab.size = 1.7
  par(mfrow=c(1,1))
  theta.stars.m = matrix(unlist(theta.stars),
                         nrow=length(theta.stars),ncol=length(theta.stars[[1]]),
                         byrow=TRUE)
  beta.paths = apply(theta.stars.m,2,get.summ.stat,getpm,lambdas)
  lgrid = beta.paths[[1]]$lgrid
  coefpaths = lapply(beta.paths,function(x){x$summ.stat})
  coefpaths.m = matrix(unlist(coefpaths),nrow=length(lgrid),ncol=dim.s$optim.exp)
  #coefpaths.m = coefpaths.m - tile(apply(coefpaths.m,2,min),dim(coefpaths.m))
  #browser()
  #lgrid = lgrid[2:length(lgrid)]
  #coefpaths.m = coefpaths.m[2:length(lgrid),]
  ylims = c(min(c(coefpaths.m,par$beh)),max(c(coefpaths.m,par$beh)))
  xlims = c(min(lgrid),max(lgrid))
  alpha=0.1
  plot(lgrid,coefpaths.m[,1],ylim=ylims,xlim=xlims,type='l',
       col=rgb(red=0,green=0,blue = 0, alpha=alpha),
       xlab = TeX("$\\lambda$"),ylab="",lty=1,cex.lab=cex.lab.size,
       panel.first = c(abline(h = par$beh,lty=1:dim(coefpaths.m)[2],
                              col=1:dim(coefpaths.m)[2]) ),log='x',cex=1)
  #legend(bquote(alpha == .(value)))#,log='y')
  legend("topright",legend=c(1:dim(coefpaths.m)[2]),
         lty=c(1:dim(coefpaths.m)[2]),col=c(1:dim(coefpaths.m)[2]))
  for(j in 2:dim(coefpaths.m)[2]){
    lines(lgrid,coefpaths.m[,j],lty=j,col=j,cex=1)#,log='x')
  }
}

lamb.plot17=function(theta.stars,parbs,lambdas){
  # bottom left top right
  par(mar=c(4,6,1,1))
  cex.lab.size = 1.7
  par(mfrow=c(1,1))
  theta.stars.m = unlist.reshape(theta.stars)
  parbs.m = unlist.reshape(parbs)
  sds = unlist(lapply(exps,function(x){x$input.ep.seed}))
  ####
  df = cbind(theta.stars.m-parbs.m)
  df[sds==1,]
  lm.u.s = sort(unique(lambdas))
  get.first = function(x,lm.u.s){
  
    for (i in 1:length(x)){
      first=i
      if (abs(x[i])<1e-3){
        break
      }
      
    }
    lm.u.s[first]
  }
  firsts=list()
  for (j in 1:length(unique(sds))){
    firsts[[j]] = apply(df[sds==j,],2,get.first,lm.u.s=lm.u.s)
  }

  fm = matrix(unlist(firsts),nrow=length(unique(sds)),ncol=dim.s$optim.exp,byrow=TRUE)
  #hist(fm[,1],xlim=c(min(fm),max(fm)),col=1,xlab=TeX("$\\lambda_{i}^*$"),main="")
  cex.lab.size = 1.7

  ylimt=3.5#max(density(fm[,2])$y,from=min(lambdas),to=max(lambdas))

  densities = list()
  for (i in 1:dim.s$optim.exp){
    densities[[i]] = density(fm[,i],from=min(lambdas),to=max(lambdas))$y
  }
  myde = unlist(densities)
  ylimt = max(myde)
  browser()
  fmd = density(fm[,2],from=min(lambdas),to=max(lambdas))
  
  xlim=c(fmd$x[2],max(lambdas))
  xlim = c(min(lambdas),max(lambdas))
  #fmd = density(fm[,1],from=min(lambdas),to=max(lambdas))

  plot(fmd$x,fmd$y,ylim=c(-1,ylimt),xlim=xlim,col=1,ylab="Density.test",
       xlab=TeX("$\\zeta_{k}^*$"),main="",type='l',cex.lab=cex.lab.size,
       cex=1,lty=2,log='x')
  

  if (dim.s$optim.exp>2){
  if (plot.17){
      fmd = density(fm[,17],from=min(lambdas),to=max(lambdas))
      #lines(fmd$x,fmd$y,lty=1,col="white",lwd=2)
      lines(fmd$x,fmd$y,lty=4,col=17,lwd=2)
  }
  for (i in 3:dim.s$optim.exp){
    #hist(fm[,i],add=TRUE,col=i,xlab="")
    fmd = density(fm[,i],from=min(lambdas),to=max(lambdas))
    if(i==17){
      if(plot.17){
        #lines(fmd$x,fmd$y,lty=2)
      }else{
        lines(fmd$x,fmd$y,lty=2)
      }
    }else{
    lines(fmd$x,fmd$y,lty=2)
    }
  }
  }
  #legend("topleft",legend=paste0("i=",1:dim.s$optim.exp),lty=1:dim.s$optim.exp,
  #       col=1:dim.s$optim.exp)
  if(plot.17){
  legend("topright",legend=c("k = 2",
                             expression(paste("k" != 2,", 17")),
                             "k = 17"),lty=c(1,2,4),
         col=c(2,1,17),lwd=c(2,1,2))
  }else{
    legend("topright",legend=c("k = 2",expression("k"!= 2)),lty=c(1,2),
           col=c(2,1),lwd=c(2,1))
    
    
  }
  fmd = density(fm[,2],from=min(lambdas),to=max(lambdas))
  lines(fmd$x,fmd$y,lty=1,col=2,lwd=2)
  #legend("topright",legend=c(expression(SSE[V[lambda]]),
  #                           expression(paste("|","!=","*",-V[lambda],"|"))),lty=1:2)
  
}



ix.info = function(ix){
  sbar=get.s.stat(exps[[ix]]$eps.raw,mean)
  abar=get.a.stat(exps[[ix]]$eps.raw,mean)
  rbar = get.r.stat(exps[[ix]]$eps.raw,mean)
  hist(unlist(lapply(exps[[ix]]$eps.v.info,function(x){x$ratio})),
       main= paste0("lambda=",lambdas[ix]),xlab='pbeta/pb')
  plot(exps[[ix]]$optim.res$par,
       main=paste0("val",round(exps[[ix]]$theta.star.value.train,4),
                   ";sbar", round(sbar,3)
),
       xlab="beta")
  plot(exps[[ix]]$mle.parb,
           main=paste0(
       # "val",round(exps[[ix]]$theta.star.value.train,4),
       # ";sbar", round(sbar,3),
       # ";abar", round(abar,3),
       # ";rbar", round(rbar,3)
       ";abar", round(abar,3),
       ";rbar", round(rbar,3)
            ),
       xlab="b")
  list(lambda=exp.list[[ix]]$lambda,
       sbar=sbar,
       abar=abar,
       rbar = rbar,
       val=exps[[ix]]$theta.star.value.train,
       var.v= exps[[ix]]$var.v
       #,
       #exps[[ix]]$optim.res$par
       #,mle=exps[[ix]]$mle.parb
  )
  
}

plot.d = function(){
  par(mfrow=c(4,3))
  #  bottom, left, top, and right.
  par(mar=c(4,2,2,2))
  print("LOW VALUE")
  ixs=which(theta.star.value.train<1e-21)
  for (ix in ixs){
    ix.info(ix)
  }
  print("HIGH VALUE")
  ixs = which(theta.star.value.train>.2)
  for (ix in ixs){
    ix.info(ix)
  }
  
}



  make.scale.plot = function(theta.star.value.train,label){
    print("need to deal with sigma")
    v.train.med = get.summ.stat(theta.star.value.train,median,lambdas)
    lgrid = v.train.med$lgrid
    v.train.med = v.train.med$summ.stat
    v.train.sse = get.summ.stat(theta.star.value.train,
                                function(x){sqrt(var(x))},lambdas)$summ.stat
    
    v.train.med = (v.train.med - mean(v.train.med))/sd(v.train.med)
    v.train.sse = (v.train.sse - mean(v.train.sse))/sd(v.train.sse)
    val.dev =  abs(max(v.train.med)-v.train.med)
    # plot(lgrid,val.dev,type='l')
    scale.plot(lgrid,v.train.sse,val.dev,label)
  }

  

scale.plot = function(ls,x,y,label,plott=TRUE,pick.lam=FALSE){
  cex.lab.size = 1.7
  #par(mfrow=c(1,1))
  #plot(ls,x,type='l',xlab='')
  #plot(ls,y,type='l',lty=2)
  
  #plot(ls,x,type='l',xlab='')
  #lines(ls,y,lty=2)
  
  scale = function(x,xs,ys){
    slope = ((min(xs)-max(xs))/(min(ys)-max(ys)))
    s.y = (x-(min(xs)-min(ys)*slope))/slope
    s.y
  }
  
  s.x = scale(x,x,y)
  if (plott){
  plot(ls,s.x,type='l',xlab=TeX("$\\lambda$"),ylab=expression(paste(sigma," or ","|",V,"*",-V[lambda],"|")),cex.lab=cex.lab.size,main=label) #,yaxt='n')
  #print(s.x)
  #print(y)
  lines(ls,y,lty=2)
  
  legend("topleft",legend=c(expression(sigma[V[lambda]]),expression(paste("|",V,"*",-V[lambda],"|"))),lty=1:2)
  #plot(ls,s.x+y,type='l',xlab=TeX("$\\lambda$"))
  }
  #lines(ls,x+y,lty=2)
  ix = which((s.x+y)==min(s.x+y))
  diff=abs(s.x-y)
  ix = which(diff==min(diff))
  
  #print(ls[ix])
  #text(ls[ix],0.1,expression(paste(lambda,"*")))
  
  if(pick.lam){
  abline(v=ls[ix],lty=3)
  mtext(expression(paste("*")),side=1,at=ls[ix])
  }
  print(c("PICKED LAMBDA=",ls[ix]))
  
  list(raw.which = which((x+y)==min(x+y)),
       sca.which=which((s.x+y)==min(s.x+y)),
       ix.star=ix)
}

