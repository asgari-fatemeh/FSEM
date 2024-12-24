#################extracting norm values##################
library(Matrix)
library(tidyverse)
library(data.table)
library(writexl)



dat_list <-results


n<-length(dat_list)
n.coef.fac<-length(dat_list[[1]]$model.fit$var$indicators)
n.coef.sem<-length(dat_list[[1]]$model.fit$var$factors)
nam.fac<-dat_list[[1]]$model.fit$var$indicators
nam.sem<-dat_list[[1]]$model.fit$var$latents
n.sigma.fac.vec<-n.coef.fac
n.sigma.fac.val<-n.coef.fac
n.sigma.sem.vec<-n.coef.sem
n.sigma.sem.val<-n.coef.sem
n.intercept<-n.coef.fac
n.sigma.err<-n.coef.fac
no.f<-length(dat_list[[1]]$model.fit$mod$factorModel)
norm<-list()
for (i in 1:n) {
  res <- dat_list[[i]]
  N<-length(unique(res$simulation$data$.id)) 
  M<-max(sapply(1:N,function(j){max(sapply(1:length(res$model.fit$mod$factorModel), function(jj){length(res$simulation$data$.t[res$simulation$data$.id==j&res$simulation$data$.ind==jj])}))}))
  n.samp<-N
  n.t<-M
  n.eigen.fac<-length(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std)
  for (k in 1:n.eigen.fac) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]])
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][kk,1]*
        res$simulation$params.org.eval$sigma.fac.eigvec.std[[k]][kk,1]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,1]<--1*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,1]
    }
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][kk,2]*
        res$simulation$params.org.eval$sigma.fac.eigvec.std[[k]][kk,2]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,2]<--1*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,2]
    }
  }
  
  n.eigen.sem<-length(res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std)
  for (k in 1:n.eigen.sem) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]])
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][kk,1]*
        res$simulation$params.org.eval$sigma.sem.eigvec.std[[k]][kk,1]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,1]<--1*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,1]
    }
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][kk,2]*
        res$simulation$params.org.eval$sigma.sem.eigvec.std[[k]][kk,2]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,2]<--1*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,2]
    }
  }
  
  norm[[paste0("res",i)]]<-list(N=N,M=M)
  for (j in 1:n.coef.fac) {
    namm<-names(res$estimation$result$params.estimated.eval$coef.fac.std[[j]])
  if(length(res$estimation$result$params.estimated.eval$coef.fac.std[[j]])>0){
    m<-length(namm)
    for (k in 1:m) {
      for (kk in 1:length(res$model.fit$mod$factorModel)) {
       if(res$model.fit$mod$factorModel[[kk]]$indicator==nam.fac[j]){
         nn<-kk
       } 
        for (kk2 in 1:length(res$model.fit$mod$factorModel[[k]]$factor)) {
          if(res$model.fit$mod$factorModel[[kk]]$factor[kk2]==namm[k]){
            mm<-kk2
          } 
        }
      }
      if(res$model.fit$mod$factorModel[[nn]]$effect[mm]=="historical"){
        norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]<-1/(200^2)*sum((res$estimation$result$params.estimated.eval$coef.fac.std[[j]][[k]]-
                                                                                res$simulation$params.org.eval$coef.fac.std[[j]][[k]])^2)
      }else{
        norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]<-1/200*sum((res$estimation$result$params.estimated.eval$coef.fac.std[[j]][[k]]-
                                                                                res$simulation$params.org.eval$coef.fac.std[[j]][[k]])^2)
      }
    }
  }  
  }
  for (j in 1:n.coef.sem) {
    namm<-names(res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    if(length(namm)>0){
      m<-length(namm)
      for (k in 1:m) {
        for (kk in 1:length(res$model.fit$mod$regression)) {
          if(res$model.fit$mod$regression[[kk]]$response==nam.sem[j]){
            nn<-kk
          } }
          for (kk2 in 1:length(res$model.fit$mod$regression[[nn]]$covariate)) {
            if(res$model.fit$mod$regression[[kk]]$covariate[kk2]==namm[k]){
              mm<-kk2
            } 
          }

        if(res$model.fit$mod$regression[[nn]]$effect[mm]=="historical"){
          norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]]<-1/(200^2)*sum((res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.sem.std[[nam.sem[j]]][[k]])^2)
        }else{
          norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]]<-1/200*sum((res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.sem.std[[nam.sem[j]]][[k]])^2)
        }
      }
    }
  }
  for (j in 1:n.intercept) {
    norm[[paste0("res",i)]]$intercept[[nam.fac[j]]]<-1/200*sum((res$estimation$result$params.estimated.eval$intercept.std[[j]]-
                                          res$simulation$params.org.eval$intercept.std[[j]])^2)
  }
  for (j in 1:n.sigma.fac.vec) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[nam.fac[j]]])[2]
    norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]][k]<-1/200*sum((res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[j]][,k]-
                                                        res$simulation$params.org.eval$sigma.fac.eigvec.std[[j]][,k])^2)  
    }
  }
  for (j in 1:n.sigma.sem.vec) {
    norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]][k]<-1/200*sum((res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[j]][,k]-
                                                        res$simulation$params.org.eval$sigma.sem.eigvec.std[[j]][,k])^2)  
    } 
  }
  for (j in 1:n.sigma.fac.val) {
    norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]][k]<-(res$estimation$result$params.estimated.eval$sigma.fac.eigval.std[[j]][k]-
                                            res$simulation$params.org.eval$sigma.fac.eigval.std[[j]][k])^2  
    }
  }
  for (j in 1:n.sigma.sem.val) {
    norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]][k]<-(res$estimation$result$params.estimated.eval$sigma.sem.eigval.std[[j]][k]-
                                            res$simulation$params.org.eval$sigma.sem.eigval.std[[j]][k])^2  
    } 
  }
  for (j in 1:n.sigma.err) {
    norm[[paste0("res",i)]]$sigma.error[[nam.fac[j]]]<-(res$estimation$result$params.estimated$totparam$sigma.error[[j]]-
      res$simulation$params.org.eval$sigma.error[[j]])^2
  }
}


#############extracting mse values###############################
comb<-expand.grid(N=c(50,100),M=c(10,20))
comb<-comb[order(comb$N),]
mse<-list()
ress<-dat_list[[1]]
for (l in 1:dim(comb)[1]) {
  mse[[l]]<-list(N=comb$N[l],M=comb$M[l])
  for (j in 1:n.coef.fac) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
    mse[[l]][[nam.fac[j]]]<-data.frame(NA)
    if(length(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])>0){
      m<-length(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
      for (k in 1:m) {
        mse[[l]][[nam.fac[j]]][[namm[k]]]<-0
        count<-0
        for (i in 1:n) {
          if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
            mse[[l]][[nam.fac[j]]][[namm[k]]]<-mse[[l]][[nam.fac[j]]][[namm[k]]]+
              norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]  
            count<-count+1
          }
        }
        mse[[l]][[nam.fac[j]]][[namm[k]]]<-1/count*mse[[l]][[nam.fac[j]]][[namm[k]]]
      }
    }
    mse[[l]][[nam.fac[j]]]$intercept<-0
    count<-0
    for (i in 1:n) {
      if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
        mse[[l]][[nam.fac[j]]]$intercept<-mse[[l]][[nam.fac[j]]]$intercept+norm[[paste0("res",i)]]$intercept[[nam.fac[j]]]
      count<-count+1
      }
    }
    mse[[l]][[nam.fac[j]]]$intercept<-1/count*mse[[l]][[nam.fac[j]]]$intercept
    for (k in 1:2) {
      mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]<-0
      mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]<-0
      count<-0
      for (i in 1:n) {
        if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
        mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]<-mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]+
          norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]][k]
        mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]<-mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]+
          norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]][k]
        count<-count+1
        }
      }
      mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]<-1/count*mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]
      mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]<-1/count*mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]
    }
    mse[[l]][[nam.fac[j]]]$sig.err<-0
    count<-0
    for (i in 1:n) {
      if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
      mse[[l]][[nam.fac[j]]]$sig.err<-mse[[l]][[nam.fac[j]]]$sig.err+norm[[paste0("res",i)]]$sigma.err[[nam.fac[j]]] 
      count<-count+1
      }
    }
    mse[[l]][[nam.fac[j]]]$sig.err<-1/count*mse[[l]][[nam.fac[j]]]$sig.err
  } 
  
for (j in 1:n.coef.sem) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    mse[[l]][[nam.sem[j]]]<-data.frame(NA)
    if(length(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])>0){
      m<-length(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
      for (k in 1:m) {
        mse[[l]][[nam.sem[j]]][[namm[k]]]<-0
        count<-0
        for (i in 1:n) {
          if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
          mse[[l]][[nam.sem[j]]][[namm[k]]]<-mse[[l]][[nam.sem[j]]][[namm[k]]]+
            norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]] 
          count<-count+1
          }
        }
        mse[[l]][[nam.sem[j]]][[namm[k]]]<-1/count*mse[[l]][[nam.sem[j]]][[namm[k]]]
      }
    }
    for (k in 1:2) {
      mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]<-0
      mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]<-0
      count<-0
      for (i in 1:n) {
        if(norm[[paste0("res",i)]]$N==comb[l,]$N&norm[[paste0("res",i)]]$M==comb[l,]$M){
        mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]<-mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]+
          norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]][k]
        mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]<-mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]+
          norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]][k]
        count<-count+1
        }
      }
      mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]<-1/count*mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]
      mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]<-1/count*mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]
    }
} 
  mse[[l]]$numb.runs<-count
}

#################extracting tables###############################
tabl.fac<-data.frame(N=comb$N,M=comb$M)
for (l in 1:dim(comb)[1]) {
  for (j in 1:n.coef.fac) {
    tabl.fac[[paste0("beta",j)]][l]<-mse[[l]][[nam.fac[j]]]$intercept
    namm<-names(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
    for (m in 1:length(namm)) {
      tabl.fac[[paste0("lambda",j,m)]][l]<-mse[[l]][[nam.fac[j]]][[namm[m]]]
    }
    tabl.fac[[paste0("phi",j,1)]][l]<-mse[[l]][[nam.fac[j]]]$eig.fuc.1
    tabl.fac[[paste0("phi",j,2)]][l]<-mse[[l]][[nam.fac[j]]]$eig.fuc.2
    tabl.fac[[paste0("nu",j,1)]][l]<-mse[[l]][[nam.fac[j]]]$eig.val.1
    tabl.fac[[paste0("nu",j,2)]][l]<-mse[[l]][[nam.fac[j]]]$eig.val.2
    tabl.fac[[paste0("sigma",j)]][l]<-mse[[l]][[nam.fac[j]]]$sig.err
  }
}


tabl.sem<-data.frame(N=comb$N,M=comb$M)
for (l in 1:dim(comb)[1]) {
  for (j in 1:n.coef.sem) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    if(length(namm)>0){
      for (m in 1:length(mse[[l]][[nam.sem[j]]])) {
        tabl.sem[[paste0("gamma",j,m)]][l]<-mse[[l]][[nam.sem[j]]][[namm[m]]]
      } 
    }
    tabl.sem[[paste0("psi",j,1)]][l]<-mse[[l]][[nam.sem[j]]]$eig.fuc.1
    tabl.sem[[paste0("psi",j,2)]][l]<-mse[[l]][[nam.sem[j]]]$eig.fuc.2
    tabl.sem[[paste0("mu",j,1)]][l]<-mse[[l]][[nam.sem[j]]]$eig.val.1
    tabl.sem[[paste0("mu",j,2)]][l]<-mse[[l]][[nam.sem[j]]]$eig.val.2
  }
} 






