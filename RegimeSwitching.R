
### Regime Switch Simu

sim.RS <- function(P, NSimu, µ, rho, sigma){
  draw<-runif(NSimu)
  state<-c()
  SIni<-sample.int(dim(P)[1], 1)
  XIni<-µ[SIni]
  XSim<-c()
  for (i in 1:NSimu){
    if (i==1){
      for (j in 1:dim(P)[1]){
        if (draw[i]<sum(P[SIni,][1:j])){
          Fstate<-j
          state<-c(state, Fstate)
          FX<-µ[state[i]]+rho[state[i]]*XIni+sigma[state[i]]*rnorm(1,0,1)
          XSim<-c(XSim, FX)
          break
        }
      }
    }
    else{
      for (j in 1:dim(P)[1]){
        if (draw[i]<sum(P[state[i-1],][1:j])){
          Fstate<-j
          state<-c(state, Fstate)
          break
        }
      }
    }
    FX<-µ[state[i]]+rho[state[i]]*XSim[i-1]+sigma[state[i]]*rnorm(1,0,1)
    XSim<-c(XSim, FX)
  }
  Fin<-list('XSim'=XSim, 'XIni'=XIni, 'state'=state, 'SIni'=SIni)
  return(Fin)
}

### Hamilton filter

param2model <- function(param,K){#Not In Use
  P <- matrix(0,K,K)
  P[,1] <- exp(param[1:K])/(1+exp(param[1:K]))
  
  for(j in 2:(K-1)){
    param_j <- param[(K*(j-1)+1):(K*j)]
    sum_previous_Proba <- matrix(P[,1:(j-1)],K,j-1) %*% matrix(1,j-1,1)
    P[,j] <- (1-sum_previous_Proba)*exp(param_j)/(1+exp(param_j))
  }
  sum_previous_Proba <- matrix(P[,1:(K-1)],K,K-1) %*% matrix(1,K-1,1)
  P[,K] <- 1-sum_previous_Proba
  return(P)
}

param2model<-function(param,K){
  mu<-param[7:9]
  phi<-abs(param[10:12])
  sigmoid<-param[13:15]
  MC<-c()
  for (i in 1:K){
    DoFProba<-abs(1-sum(abs(param[(i*(K-1)-(K-2)):(i*(K-1))])))
    rowi<-(c(abs(param[(i*(K-1)-(K-2)):(i*(K-1))]),DoFProba))/sum(c(abs(param[(i*(K-1)-(K-2)):(i*(K-1))]),DoFProba))
    MC<-rbind(MC,rowi)
    if (phi[i]>0.99){
      phi[i]<-0.99
    }
    sigmoid[i]<-max(sigmoid[i],(1/(sqrt(2*pi))+0.01))
  }
  fin<-list('MC'=MC,'mu'=mu,'phi'=phi,'sigmoid'=sigmoid)
  return(fin)
}

paramToModel<-function(param,K){#Not In Use
  param<-abs(param)
  for (i in 1:K){
    if (i==1){
      i2K<-1-sum(param[1:(K-1)])
      MC<-c(param[1:(K-1)],i2K)
    }
    else{
      i2K<-1-sum(param[((i-1)*(K-1)+1):(i*(K-1))])
      rowi<-c(param[((i-1)*(K-1)+1):(i*(K-1))],i2K)
      MC<-rbind(MC,rowi)
    }
  }
  return(MC)
}

ParamToModelNODOF<-function(param,K){#NOT IN USE
  param<-abs(param)
  for (i in 1:K){
    if (i==1){
      row1<-param[1:K]/(sum(param[1:K]))
      MC<-c(row1)
    }
    else{
      rowi<-param[((i-1)*K+1):(i*K)]/(sum(param[((i-1)*K+1):(i*K)]))
      MC<-rbind(MC,rowi)
    }
  }
  return(MC)
}

## Optim 3 regimes

RSLik3<-function(params, x){
  P<-param2model(params,3)$MC
  mu<-param2model(params,3)$mu
  phi<-param2model(params,3)$phi
  sigmoid<-param2model(params,3)$sigmoid
  eta<-c()
  KsiIni<-(P%^%100)[1,]
  for (i in 1:length(x)){
    etai<-c()
    for (j in 1:dim(P)[1]){
      if(i==1){
        etaij<-(1/(sigmoid[j]*(2*pi)**(1/2)))*exp(((-1)/2)*((x[i]-(mu[j]/(1-phi[j])))/(sigmoid[j]))**2)
        etai<-c(etai, etaij)
      }
      else{
        etaij<-(1/(sigmoid[j]*(2*pi)**(1/2)))*exp(((-1)/2)*((x[i]-(mu[j]+phi[j]*x[i-1]))/(sigmoid[j]))**2)
        etai<-c(etai, etaij)
      }
    }
    if (i==1){
      eta<-c(eta, etai)
      ksi1<-((t(P)%*%KsiIni)*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%KsiIni)*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%KsiIni)*etai)))
      Lik<-c(Liki)
      ksi<-c(t(ksi1))
    }
    else if (i==2){
      eta<-rbind(eta, etai)
      ksi2<-((t(P)%*%ksi)*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi)*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi)*etai)))
      Lik<-c(Lik, Liki)
      ksi<-rbind(ksi, t(ksi2))
    }
    else{
      eta<-rbind(eta, etai)
      ksit<-((t(P)%*%ksi[i-1,])*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi[i-1,])*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi[i-1,])*etai)))
      Lik<-c(Lik, Liki)
      ksi<-rbind(ksi, t(ksit))
    }
  }
  AllLik<-sum(log(Lik))
  return(-AllLik)
}

##Lik Val

RSLik3Val<-function(params, x){
  P<-param2model(params,3)$MC
  mu<-param2model(params,3)$mu
  phi<-param2model(params,3)$phi
  sigmoid<-param2model(params,3)$sigmoid
  eta<-c()
  KsiIni<-(P%^%100)[1,]
  for (i in 1:length(x)){
    etai<-c()
    for (j in 1:dim(P)[1]){
      if(i==1){
        etaij<-(1/(sigmoid[j]*(2*pi)**(1/2)))*exp(((-1)/2)*((x[i]-(mu[j]/(1-phi[j])))/(sigmoid[j]))**2)
        etai<-c(etai, etaij)
      }
      else{
        etaij<-(1/(sigmoid[j]*(2*pi)**(1/2)))*exp(((-1)/2)*((x[i]-(mu[j]+phi[j]*x[i-1]))/(sigmoid[j]))**2)
        etai<-c(etai, etaij)
      }
    }
    if (i==1){
      eta<-c(eta, etai)
      ksi1<-((t(P)%*%KsiIni)*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%KsiIni)*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%KsiIni)*etai)))
      Lik<-c(Liki)
      ksi<-c(t(ksi1))
    }
    else if (i==2){
      eta<-rbind(eta, etai)
      ksi2<-((t(P)%*%ksi)*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi)*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi)*etai)))
      Lik<-c(Lik, Liki)
      ksi<-rbind(ksi, t(ksi2))
    }
    else{
      eta<-rbind(eta, etai)
      ksit<-((t(P)%*%ksi[i-1,])*etai)/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi[i-1,])*etai))))
      Liki<-c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi[i-1,])*etai)))
      Lik<-c(Lik, Liki)
      ksi<-rbind(ksi, t(ksit))
    }
  }
  AllLik<-sum(log(Lik))
  Fin<-list('Lik'=-AllLik, 'eta'=eta, 'ksi'=ksi, 'KsiIni'=KsiIni)
  return(Fin)
}

## Optim function

MonteCarloMLRS<-function(Obs,iter){
  Val<-c()
  for (i in 1:(length(Obs))){
    if (Obs[i]<0){
      Setting<-'Neg'
      break
    }
    else{
      Setting<-'Pos'
    }
  }
  if (Setting=='Neg'){
    for (i in 1:iter){
      TempParams<-c(runif(6,0,1),mean(Obs),-mean(Obs),runif(1,min(mean(Obs),-mean(Obs)),max(mean(Obs),-mean(Obs))),runif(3,0,0.99),runif(3,1,10))
      TempLik<-RSLik3(TempParams,Obs)
      if (is.nan(TempLik)==FALSE){
        TempParams<-c(TempParams,TempLik)
        Val<-rbind(Val,TempParams)
      }
    }
  }
  else{
    for (i in 1:iter){
      TempParams<-c(runif(6,0,1),runif(3,0,mean(Obs)),runif(3,0,0.99),runif(3,1,10))
      TempLik<-RSLik3(TempParams,Obs)
      if (is.nan(TempLik)==FALSE){
        TempParams<-c(TempParams,TempLik)
        Val<-rbind(Val,TempParams)
      }
    }
  }
  return(Val)
}

EstimRS3AR1<-function(Obs,iter){
  a<-MonteCarloMLRS(Obs,iter)
  MCParams<-a[FindIndex(a[,(dim(a)[2])],min(a[,(dim(a)[2])])),][1:15]
  OptimPar<-optim(MCParams,
                  RSLik3,
                  x=Obs,
                  control=list(trace=FALSE,maxit=10000))
  return(OptimPar$par)
}

## Kim filter

CompKimKsi<-function(ksis, P){
  for (i in 1:(dim(ksis)[1])){
    if (i==1){
      KimKsi<-ksis[dim(ksis)[1],]
    }
    else if (i==2){
      KimKsit<-c(ksis[((dim(ksis)[1])-(i-1)),]*(P%*%(KimKsi/(t(P)%*%ksis[((dim(ksis)[1])-(i-1)),]))))
      KimKsi<-rbind(KimKsit, KimKsi)
    }
    else{
      KimKsit<-c(ksis[((dim(ksis)[1])-(i-1)),]*(P%*%(KimKsi[1,]/(t(P)%*%ksis[((dim(ksis)[1])-(i-1)),]))))
      KimKsi<-rbind(KimKsit, KimKsi)
    }
  }
  return(KimKsi)
}

### Forecast

##Forecast by multiple path simulation 3 regimes

sim.RSFcast <- function(P, NSimu, µ, rho, sigma, SIni,XIni){
  draw<-runif(NSimu)
  state<-c()
  XSim<-c()
  for (i in 1:NSimu){
    if (i==1){
      for (j in 1:dim(P)[1]){
        if (draw[i]<sum(P[SIni,][1:j])){
          Fstate<-j
          state<-c(state, Fstate)
          FX<-µ[state[i]]+rho[state[i]]*XIni+sigma[state[i]]*rnorm(1,0,1)
          XSim<-c(XSim, FX)
          break
        }
      }
    }
    else{
      for (j in 1:dim(P)[1]){
        if (draw[i]<sum(P[state[i-1],][1:j])){
          Fstate<-j
          state<-c(state, Fstate)
          break
        }
      }
    }
    FX<-µ[state[i]]+rho[state[i]]*XSim[i-1]+sigma[state[i]]*rnorm(1,0,1)
    XSim<-c(XSim, FX)
  }
  Fin<-list('XSim'=XSim, 'XIni'=XIni, 'state'=state, 'SIni'=SIni)
  return(Fin)
}

ForecastInt<-function(params,x,steps,confidence,iterations){
  path<-c()
  upbound<-c()
  lowbound<-c()
  PreviousVal<-RSLik3Val(params,x)
  LikLastState<-FindIndex(PreviousVal$ksi[dim(PreviousVal$ksi)[1],],max(PreviousVal$ksi[dim(PreviousVal$ksi)[1],]))
  for (i in 1:iterations){
    pathi<-sim.RSFcast(param2model(params,3)$MC,steps,param2model(params,3)$mu,param2model(params,3)$phi,param2model(params,3)$sigmoid,LikLastState,x[length(x)])
    path<-rbind(path,pathi$XSim)
  }
  for (i in 1:steps){
    upboundi<-quantile(path[,i],(1-(confidence/2)))
    upbound<-c(upbound,upboundi)
    lowboundi<-quantile(path[,i],confidence/2)
    lowbound<-c(lowbound,lowboundi)
  }
  Fin<-list("UpBound"=upbound, "LowBound"=lowbound)
  return(Fin)
}


##Point Forecast

ForecastPointRS<-function(params,x,steps,iterations){
  path<-c()
  pmean<-c()
  PreviousVal<-RSLik3Val(params,x)
  LikLastState<-FindIndex(PreviousVal$ksi[dim(PreviousVal$ksi)[1],],max(PreviousVal$ksi[dim(PreviousVal$ksi)[1],]))
  for (i in 1:iterations){
    pathi<-sim.RSFcast(param2model(params,3)$MC,steps,param2model(params,3)$mu,param2model(params,3)$phi,param2model(params,3)$sigmoid,LikLastState,x[length(x)])
    path<-rbind(path,pathi$XSim)
  }
  for (i in 1:steps){
    pmeani<-mean(path[,i])
    pmean<-c(pmean,pmeani)
  }
  return(pmean)
}


## States vs ksi

TransfState<-function(data,state){
  BinState<-c()
  for (i in 1:length(data)){
    if (data[i]==state){
      BinStatei<-1
      BinState<-c(BinState, BinStatei)
    }
    else{
      BinStatei<-0
      BinState<-c(BinState, BinStatei)
    }
  }
  return(BinState)
}

###Analytical forecast

ForecastPointRSProb<-function(params,horizon,x){
  regimes<-3
  P<-param2model(params,3)$MC
  mu<-param2model(params,3)$mu
  phi<-param2model(params,3)$phi
  sigmoid<-param2model(params,3)$sigmoid
  etas<-RSLik3Val(params,x)$eta
  etai<-c()
  path<-c()
  ksi<-RSLik3Val(params,x)$ksi
  upto<-dim(ksi)[1]
  for (i in 1:horizon){
    etai<-c()
    pathij<-c()
    for (j in 1:regimes){
      subpath<-mu[j]+phi[j]*x[upto+i-1]
      lb<-ksi[upto+i-1,]%*%P[,j]
      pathij<-c(pathij, subpath*lb)
    }
    pathi<-sum(pathij)
    path<-c(path,pathi)
    x<-c(x,pathi)
    for (j in 1:regimes){
      etaij<-(1/(sigmoid[j]*(2*pi)**(1/2)))*exp(((-1)/2)*((x[upto+i]-(mu[j]+phi[j]*x[upto+i-1]))/(sigmoid[j]))**2)
      etai<-c(etai, etaij)
    }
    etas<-rbind(etas,etai)
    ksit<-((t(P)%*%ksi[upto+i-1,])*etas[upto+i,])/(c((t(rep(1, dim(P)[1]))%*%((t(P)%*%ksi[upto+i-1,])*etas[upto+i,]))))
    ksi<-rbind(ksi, t(ksit))
  }
  return(path)
}
