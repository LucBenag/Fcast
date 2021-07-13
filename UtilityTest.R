### Utilities

## FindIndex

FindIndex<-function(Vector, Value){
  Index<-c()
  for (i in 1:length(Vector)){
    if (Value==Vector[i]){
      Index<-c(Index, i)
    }
  }
  return(Index)
}

## Growth

GRate<-function(x){
  for (i in 1:length(x)){
    if (i==1){
      g<-0
    }
    else{
      gi<-((x[i]-x[i-1])/x[i-1])*100
      g<-c(g,gi)
    }
  }
  g[1]<-mean(g[-1])
  return(g)
}

## Growth to index

GrowthToValue<-function(LastVal,Growth){
  Growth<-Growth/100
  for (i in 1:length(Growth)){
    if (i==1){
      val<-LastVal*(1+Growth[i])
    }
    else{
      vali<-val[i-1]*(1+Growth[i])
      val<-c(val,vali)
    }
  }
  return(val)
}

## Percent change with corresponding previous year's month 

PYMChange<-function(index,x){
  return(((x[index]-x[index-12])/x[index-12])*100)
}

##Linear trend

DeTrend<-function(x){
  for (i in 1:length(x)){
    if (i==1){
      k<-x[1]
    }
    else{
      ki<-x[i]-x[i-1]
      k<-c(k,ki)
    }
  }
  k[1]<-mean(k)
  return(k)
}

## Num digits

IntDigits<-function(number){
  decim<-c()
  for (i in 0:10){
    deci<-10**i
    decim<-c(decim,deci)
  }
  for (i in 1:length(decim)){
    if (number>=decim[i]){
      next
    }
    return(i-1)
    break
  }
}

## PYM to val

PYMtoPointVal<-function(PYMFcast,PastVal){
  PYMFcast<-PYMFcast/100
  path<-PastVal
  for (i in 1:length(PYMFcast)){
    pathi<-path[(length(path)-11)]*(1+PYMFcast[i])
    path<-c(path,pathi)
  }
  return(path[-(1:length(PastVal))])
}

PYMtoIntVal<-function(PYMFcastUpBound,PYMFcastLowBound,PastVal){
  PYMFcastUpBound<-PYMFcastUpBound/100
  PYMFcastLowBound<-PYMFcastLowBound/100
  uppath<-PastVal
  lowpath<-PastVal
  for (i in 1:length(PYMFcastUpBound)){
    uppathi<-uppath[(length(uppath)-11)]*(1+PYMFcastUpBound[i])
    uppath<-c(uppath,uppathi)
    lowpathi<-lowpath[(length(lowpath)-11)]*(1+PYMFcastLowBound[i])
    lowpath<-c(lowpath,lowpathi)
  }
  fin<-list('UpBound'=uppath[-(1:length(PastVal))],'LowBound'=lowpath[-(1:length(PastVal))])
  return(fin)
}

### Test

##Diebold

FCastSqError<-function(PointEstim,obs){
  sqerror<-c()
  for (i in 1:length(PointEstim)){
    sqerrori<-(PointEstim[i]-obs[i])**2
    sqerror<-c(sqerror,sqerrori)
  }
  return(sqerror)
}

DieboldSeq<-function(ErrorA,ErrorB){
  Seq<-c()
  for (i in 1:length(ErrorA)){
    Seqi<-ErrorA[i]-ErrorB[i]
    Seq<-c(Seq,Seqi)
  }
  return(Seq)
}

autocov <- function(X,n){
  Steps<-length(X)
  dev1<-X[1:(Steps-n)]-mean(X)
  dev2<-X[(n+1):Steps]-mean(X)
  return((1/Steps)*sum(dev1*dev2))
}

NWAutoCov<-function(X){
  acovSeq<-c()
  for (i in 0:(length(X)-1)){
    acovi<-autocov(X,i)
    acovSeq<-c(acovSeq,acovi)
  }
  for(i in 0:(length(acovSeq)-1)){
    if (i==0){
      NeweyWest<-acovSeq[1]
    }
    else{
      NeweyWesti<-2*(1-(i/(length(acovSeq))))*acovSeq[i+1]
      NeweyWest<-c(NeweyWest,NeweyWesti)
    }
  }
  return(sum(NeweyWest))
}

DieboldTest<-function(ErrorA,ErrorB,confidence){
  acovSeq<-c()
  Seq<-DieboldSeq(ErrorA, ErrorB)
  for (i in 0:(length(Seq)-1)){
    acovi<-autocov(Seq,i)
    acovSeq<-c(acovSeq,acovi)
  }
  for(i in 0:(length(acovSeq)-1)){
    if (i==0){
      NeweyWest<-acovSeq[1]
    }
    else{
      NeweyWesti<-2*(1-(i/(length(acovSeq))))*acovSeq[i+1]
      NeweyWest<-c(NeweyWest,NeweyWesti)
    }
  }
  teststat<-sqrt(length(Seq))*((mean(Seq))/sqrt(sum(NeweyWest)))
  return(teststat)
}

## Box Jenkins NW

SimGARCHAvgNW<-function(lags,StaParams,Obs){
  SimNW<-c()
  for (i in 1:500){
    Simi<-SimGarch(lags,length(Obs),StaParams$mu,StaParams$phi,StaParams$ksi,StaParams$alpha,StaParams$beta)
    NWi<-NWAutoCov(Simi)
    SimNW<-c(SimNW,NWi)
  }
  return(mean(SimNW))
}

SimRSAvgNW<-function(StaParams,Obs){
  SimNW<-c()
  for (i in 1:500){
    Simi<-sim.RS(StaParams$MC,length(Obs),StaParams$mu,StaParams$phi,StaParams$sigmoid)$XSim
    NWi<-NWAutoCov(Simi)
    SimNW<-c(SimNW,NWi)
  }
  return(mean(SimNW))
}

## Testing interval aptitude

IntAbility<-function(UpBound,LowBound,confidence,x){
  countup<-0
  countlow<-0
  for (i in 1:length(UpBound)){
    if (UpBound[i]<x[i]){
      countup<-countup+1
    }
  }
  for (i in 1:length(LowBound)){
    if (LowBound[i]>x[i]){
      countlow<-countlow+1
    }
  }
  PropOutUp<-countup/length(UpBound)
  PropOutLow<-countlow/length(LowBound)
  cat(' Observations outside of upper bound', PropOutUp*100,'% of the time instead of', (confidence/2)*100, '% expected \n',
      'Observations outside of lower bound', PropOutLow*100,'% of the time instead of', (confidence/2)*100, '% expected \n')
  fin<-list('PropOutUp'=PropOutUp, 'PropOutLow'=PropOutLow)
  return(fin)
}

### Random Walk

SimRW<-function(sigmoid, vini, steps){
  path<-vini
  for (i in 2:steps){
    pathi<-path[i-1]+sigmoid*rnorm(1,0,1)
    path<-c(path,pathi)
  }
  return(path)
}

RWLik<-function(param, x){
  sigmoid<-abs(param[1])
  for (i in 1:length(x)){
    if (i==1){
      lik<-(1/(sigmoid*(2*pi)**(1/2)))*exp((-1/2)*(((x[i]-mean(x))/sigmoid)**2))
    }
    else{
      liki<-(1/(sigmoid*(2*pi)**(1/2)))*exp((-1/2)*(((x[i]-x[i-1])/sigmoid)**2))
      lik<-c(lik, liki)
    }
  }
  fin<-sum(log(lik))
  return(-fin)
}

EstimRW<-function(IniParam,Obs){
  OptimPar<-optim(IniParam,
                  RWLik,
                  x=Obs,
                  method='Brent',
                  lower=0,
                  upper=var(Obs))
  return(OptimPar$par)
}

FCastPointRW<-function(sigmoid,horizon,lastval){
  path<-rep(lastval,horizon)
  return(path)
}

FCastIntRW<-function(sigmoid,horizon,confidence,iterations,lastval){
  path<-c()
  upbound<-c()
  lowbound<-c()
  for (i in 1:iterations){
    pathi<-SimRW(sigmoid,lastval,horizon)
    path<-rbind(path,pathi)
  }
  for (i in 1:horizon){
    upboundi<-quantile(path[,i],(1-(confidence/2)))
    upbound<-c(upbound,upboundi)
    lowboundi<-quantile(path[,i],confidence/2)
    lowbound<-c(lowbound,lowboundi)
  }
  Fin<-list("UpBound"=upbound, "LowBound"=lowbound)
  return(Fin)
}

### Kernel Density

kernbis<-function(mu,sigma,dom){
  p<-c()
  for (i in dom){
    pri<-(1/(sigma*(pi*2)**(1/2)))*exp((-1/2)*((i-mu)/sigma)**2)
    p<-c(p,pri)
  }
  return(p)
}

KernDens<-function(x){
  kernel<-c()
  avgval<-sum(x)/length(x)
  sigval<-sqrt((sum((x-avgval)**2)/length(x)))
  domain<-(round(avgval-10*sigval):round(avgval+10*sigval))
  for (i in x){
    kerni<-kernbis(i,sigval,domain)/length(x)
    kernel<-cbind(kernel,kerni)
  }
  return(kernel)
}
