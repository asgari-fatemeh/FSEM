library(Matrix)
library(fda)

dat_list<-results
n<-length(dat_list)
mont<-100
cov.fac<-list()
lam.fac<-list()
per<-list()
fregion.bands<-list()
boot<-200
ind<-seq(1,200,by=2)
no.it<-10
cov.sem<-list()
gam.sem<-list()
per.sem<-list()
lam.fac<-list()
cov.fac<-list()
per.fac<-list()
fregion.bands.sem<-list()
fregion.bands.fac<-list()
re<-dat_list[[1]]
no.f<-length(re$model.fit$mod$factorModel)
no.r<-length(re$model.fit$mod$regression)
lat2<-re$model.fit$var$latents
n.b<-length(re$estimation$result$params.estimated$intercept[[1]])
basis<-create.bspline.basis(nbasis=n.b,rangeval = c(0,1))
times<-(seq(0,1,length.out=200))[ind]
eval1<-t(eval.basis(times,basis))
for (n.out in 1:n) {
  re<-dat_list[[n.out]]
  lam.fac[[n.out]]<-list()
  cov.fac[[n.out]]<-list()
  per.fac[[n.out]]<-list()
  fregion.bands.fac[[n.out]]<-list()
  cov.sem[[n.out]]<-list()
  gam.sem[[n.out]]<-list()
  per.sem[[n.out]]<-list()
  fregion.bands.sem[[n.out]]<-list()
  gamma1<-list()
  lambda1<-list()
  sig.sem<-list()
  samp<-unique(re$simulation$data$.id)
  n.samp<-length(samp)
  N<-n.samp
  M<-max(sapply(1:N,function(j){max(sapply(1:length(re$model.fit$mod$factorModel), function(jj){length(re$simulation$data$.t[re$simulation$data$.id==j&re$simulation$data$.ind==jj])}))}))
  n.t<-M
  for (k in 1:boot) {
    gamma1[[k]]<-list()
    lambda1[[k]]<-list()
    sig.sem[[k]]<-list()
    samples<-sort(sample(samp,replace=TRUE))
    for (j in 1:no.r) {
      lat.count<-which(lat2==re$model.fit$mod$regression[[j]]$response)
      gamma1[[k]][[lat.count]]<-list()
      if(!is.null(re$model.fit$mod$regression[[j]]$covariate)){
        la<-length(re$model.fit$mod$regression[[j]]$covariate)
        g<-list()  
        gg<-list()
        ett<-list()
        for (m in 1:mont) {
          g[[m]]<-list()
          gg[[m]]<-0
          ett[[m]]<-0
          for (i1 in samples) {
            g[[m]][[i1]]<-re$estimation$hist$hist.sem[[lat.count]]$g.eta[[m]][[i1]]
            gg[[m]]<-rbind(gg[[m]],g[[m]][[i1]])
            ett[[m]]<-c(ett[[m]],re$estimation$hist$hist.fac[[1]]$eta.l[[m]][[i1]][[lat.count]])
          }
          gg[[m]]<-gg[[m]][-1,]
          ett[[m]]<-ett[[m]][-1]
        }
        sig<-re$estimation$result$params.estimated$sigma.sem[[lat.count]]
        deltaa<-re$estimation$hist$hist.sem[[lat.count]]$deltaa
        pen.sem<-re$estimation$hist$hist.sem[[lat.count]]$pen.sem
        for (it in 1:no.it) {
          sol<-solve(kronecker(diag(samples),sig))
          gam.sum<-0
          gam1.sum<-0
          for (m in 1:mont) {
            gam<-t(gg[[m]])%*%sol%*%gg[[m]]+deltaa*pen.sem
            gam1<-t(gg[[m]])%*%sol%*%ett[[m]]
            gam.sum<-gam.sum+1/mont*gam
            gam1.sum<-gam1.sum+1/mont*gam1
          }
          gamma2.1<-solve(gam.sum)%*%gam1.sum
          sigma1<-0
          for (i1 in samples) {
            ss2<-0
            for (m in 1:mont) {
              ss<-re$estimation$hist$hist.fac[[1]]$eta.l[[m]][[i1]][[lat.count]]- re$estimation$hist$hist.sem[[lat.count]]$g[[m]][[i1]]%*%gamma2.1
              ss2<-ss2+1/mont*ss%*%t(ss)
            }
            sigma1<-sigma1+1/length(samples)*ss2
          }
          sig<-sigma1
        }
        sig.sem[[k]][[lat.count]]<-sig
        count<-0
        for (jj in 1:la) {
          count3<-which(lat2==re$model.sim$mod$regression[[j]]$covariate[[jj]])
          d.sem2<-t(eval1)%*%sig.sem[[k]][[count3]]%*%eval1 
          d.sem1<-t(eval1)%*%sig.sem[[k]][[lat.count]]%*%eval1  
          if(re$model.sim$mod$regression[[j]]$effect[[jj]]=="concurrent"){
            ev<-eval1
            ga<-gamma2.1[(count+1):(count+n.b)]
            count<-count+n.b
            gamma1[[k]][[lat.count]][[count3]]<-diag(diag(d.sem1)^(-1/2))%*%diag(diag(d.sem2)^(1/2))%*%t(ev)%*%ga
          }else{
            evv_ind<-t(eval.basis((seq(0,times[5],length.out=200))[ind],basis))
            ev<-diag(d.sem1)^(-1/2)*diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
            evv<-(diag(d.sem1)^(-1/2)*diag(d.sem2)^(1/2))[5]*kronecker(ev[,5],evv_ind)
            ga<-gamma2.1[(count+1):(count+n.b*n.b)]
            count<-count+n.b*n.b
            gamma1[[k]][[lat.count]][[count3]]<-t(evv)%*%ga
          }
        }
      }else{
        sigma1<-0
        for (i1 in samples) {
          ss2<-0
          for (m in 1:mont) {
            ss<-re$estimation$hist$hist.fac[[1]]$eta.l[[m]][[i1]][[lat.count]]  
            ss2<-ss2+1/mont*ss%*%t(ss)
          }
          sigma1<-sigma1+1/length(samples)*ss2
        }
        sig<-sigma1
        sig.sem[[k]][[lat.count]]<-sig
      }
    }
    
    for (j in 1:no.f) {
      count3<-which(lat2==re$model.sim$mod$factorModel[[j]]$factor)
      ff<-list()
      zz<-list()
      ff2<-list()
      for (m in 1:mont) {
        ff[[m]]<-0
        zz[[m]]<-0
        ff2[[m]]<-0
        for (i1 in samples) {
          ff[[m]]<-rbind(ff[[m]],re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]])
          zz[[m]]<-c(zz[[m]],re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])
          ff2[[m]]<-c(ff2[[m]],re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]])
        }
        ff[[m]]<-ff[[m]][-1,]
        zz[[m]]<-zz[[m]][-1]
        ff2[[m]]<-ff2[[m]][-1]
      }
      
      sig<-re$estimation$result$params.estimated$sigma.fac[[j]] 
      deltaa<-re$estimation$hist$hist.fac[[j]]$deltaa
      pen.fac<-re$estimation$hist$hist.fac[[j]]$pen.fac
      for (it in 1:no.it) {
        sol<-solve(kronecker(diag(samples),sig))
        lam.sum<-0
        lam1.sum<-0
        for (m in 1:mont) {
          lam<-t(ff[[m]])%*%sol%*%ff[[m]]+deltaa*pen.fac
          lam1<-t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
          lam.sum<-lam.sum+1/mont*lam
          lam1.sum<-lam1.sum+1/mont*lam1
        }
        lambda2.1<-solve(lam.sum)%*%lam1.sum
        sigma1<-0
        for (i1 in samples) {
          ss2<-0
          for (m in 1:mont) {
            ss<-as.vector(re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])-
              re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]]-re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]]%*%lambda2.1
            ss2<-ss2+1/mont*(ss%*%t(ss))
          }
          sigma1<-sigma1+1/length(samples)*ss2
        }
        sig<-sigma1
      }
      d.sem2<-t(eval1)%*%sig.sem[[k]][[count3]]%*%eval1
      if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){
        lambda1[[k]][[j]]<-as.vector(diag(d.sem2)^(1/2))
      }else{
        if(re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
          ev<-diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
          lambda1[[k]][[j]]<-t(ev)%*%lambda2.1[(n.b+1):length(lambda2.1)]
        }else{
          evv_ind<-t(eval.basis((seq(0,times[5],length.out=200))[ind],basis))
          ev<-diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
          evv<-(diag(d.sem2)^(1/2))[5]*kronecker(ev[,5],evv_ind)
          lambda1[[k]][[j]]<-t(evv)%*%lambda2.1[(n.b+1):length(lambda2.1)]
        }
      }
    }
  }
  
  for (j in 1:no.r) {
    if(!is.null(re$model.fit$mod$regression[[j]]$covariate)){
      lat.count<-which(lat2==re$model.fit$mod$regression[[j]]$response)
      cov.sem[[n.out]][[lat.count]]<-list()
      gam.sem[[n.out]][[lat.count]]<-list()
      per.sem[[n.out]][[lat.count]]<-list()
      fregion.bands.sem[[n.out]][[lat.count]]<-list()
      la<-length(re$model.fit$mod$regression[[j]]$covariate)
      count<-0
      for (jj in 1:la) {
        count3<-which(lat2==re$model.sim$mod$regression[[j]]$covariate[[jj]])
        gam.sem[[n.out]][[lat.count]][[count3]]<-list()
        cov.sem[[n.out]][[lat.count]][[count3]]<-list()
        per.sem[[n.out]][[lat.count]][[count3]]<-list()
        res<-lat2[lat.count]
        cov<-lat2[count3]
        gamma1.1<-t(sapply(1:boot,function(k){gamma1[[k]][[lat.count]][[count3]]}))
        if(re$model.sim$mod$regression[[j]]$effect[[jj]]=="concurrent"){
          gam.sem[[n.out]][[lat.count]][[count3]]<-re$estimation$result$params.estimated$gamma.param.eta.std[[lat.count]][[count3]][ind]  
          cov.sem[[n.out]][[lat.count]][[count3]]<-cov(gamma1.1)
          gam.org<-re$simulation$params.org.eval$coef.sem.std[[res]][[cov]][ind]
        }else{
          gam.sem[[n.out]][[lat.count]][[count3]]<-re$estimation$result$params.estimated$gamma.param.eta.std[[lat.count]][[count3]][ind,ind][,5]
          cov.sem[[n.out]][[lat.count]][[count3]]<-cov(gamma1.1)
          gam.org<-re$simulation$params.org.eval$coef.sem.std[[res]][[cov]][ind,ind][,5]
        }
        gam.sem[[n.out]][[lat.count]][[count3]]<-as.vector(gam.sem[[n.out]][[lat.count]][[count3]])
        cov.sem[[n.out]][[lat.count]][[count3]]<-as.matrix(cov.sem[[n.out]][[lat.count]][[count3]])
        fregion.bands.sem[[n.out]][[lat.count]][[count3]] = fregion::fregion.band(gam.sem[[n.out]][[lat.count]][[count3]],cov.sem[[n.out]][[lat.count]][[count3]], type = c("BEc"),conf.level=c(0.95))
        per.sem[[n.out]][[lat.count]][[count3]]<-mean(fregion.bands.sem[[n.out]][[lat.count]][[count3]][,3]<gam.org&gam.org<fregion.bands.sem[[n.out]][[lat.count]][[count3]][,2])
      } 
    }
  }
  for (j in 1:no.f) {
    lambda1.1<-as.matrix(t(sapply(1:boot,function(k){lambda1[[k]][[j]]})))
    count3<-which(lat2==re$model.fit$mod$factorModel[[j]]$factor)
    res<-re$model.fit$var$indicators[j]
    cov<-lat2[count3]
    if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){
      lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind]
      cov.fac[[n.out]][[j]]<-cov(lambda1.1)
      lam.org<-re$simulation$params.org.eval$coef.fac.std[[res]][[cov]][ind]
    }else{
      if(re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
        lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind]  
        cov.fac[[n.out]][[j]]<-cov(lambda1.1)
        lam.org<-re$simulation$params.org.eval$coef.fac.std[[res]][[cov]][ind]
      }else{
        lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind,ind][,5]
        cov.fac[[n.out]][[j]]<-cov(lambda1.1)
        lam.org<-re$simulation$params.org.eval$coef.fac.std[[res]][[cov]][ind,ind][,5]
      } 
    }
    lam.fac[[n.out]][[j]]<-as.vector(lam.fac[[n.out]][[j]])
    cov.fac[[n.out]][[j]]<-as.matrix(cov.fac[[n.out]][[j]])
    fregion.bands.fac[[n.out]][[j]] = fregion::fregion.band(lam.fac[[n.out]][[j]],cov.fac[[n.out]][[j]],type = c("BEc"),conf.level=c(0.95))
    per.fac[[n.out]][[j]]<-mean(fregion.bands.fac[[n.out]][[j]][,3]<lam.org&lam.org<fregion.bands.fac[[n.out]][[j]][,2])  
  }     
}


comb<-expand.grid(N=c(50,100),M=c(10,20))
comb<-comb[order(comb$N),]
cov.table.fac<-data.frame(N=comb$N,M=comb$M)

for (i in 1:dim(cov.table.fac)[1]) {
  for (j in 1:no.f) {
    cov.table.fac[[paste0("lambda",j)]][i]<-0
    count<-0
    for (k in 1:n) {
      re<-dat_list[[k]]
      N<-length(unique(re$simulation$data$.id)) 
      M<-max(sapply(1:N,function(j){max(sapply(1:length(re$model.fit$mod$factorModel), function(jj){length(re$simulation$data$.t[re$simulation$data$.id==j&re$simulation$data$.ind==jj])}))}))
      if(N==comb[i,]$N&M==comb[i,]$M){
        cov.table.fac[[paste0("lambda",j)]][i]<-cov.table.fac[[paste0("lambda",j)]][i]+
          per.fac[[k]][[j]]
        count<-count+1
      }
    } 
    cov.table.fac[[paste0("lambda",j)]][i]<-1/count*cov.table.fac[[paste0("lambda",j)]][i]
  }
}

cov.table.sem<-data.frame(N=comb$N,M=comb$M)
for (i in 1:dim(cov.table.sem)[1]) {
  for (j in 1:no.r) {
    lat.count<-which(lat2==re$model.fit$mod$regression[[j]]$response)
    la<-length(re$model.sim$mod$regression[[j]]$covariate)
    if(la>0)
      for (jj in 1:la) {
        count3<-which(lat2==re$model.sim$mod$regression[[j]]$covariate[[jj]])
        cov.table.sem[[paste0("gamma",lat.count,count3)]][i]<-0
        count<-0
        for (k in 1:n) {
          re<-dat_list[[k]]
          ree<-dat_list[[k]]
          N<-length(unique(re$simmulation$data$.id)) 
          M<-max(sapply(1:N,function(j){max(sapply(1:length(re$model.fit$mod$factorModel), function(jj){length(re$simulation$data$.t[re$simulation$data$.id==j&re$simulation$data$.ind==jj])}))}))
          if(N==comb[i,]$N&M==comb[i,]$M){
            if(length(count3)>0){
              if(!is.null(per.sem[[lat.count]][[count3]])){
                cov.table.sem[[paste0("gamma",lat.count,count3)]][i]<-cov.table.sem[[paste0("gamma",lat.count,count3)]][i]+
                  per.sem[[lat.count]][[count3]] 
              }else{
                cov.table.sem[[paste0("gamma",lat.count,count3)]][i]<-NA
              } 
              count<-count+1
            }
          }
        }
        cov.table.sem[[paste0("gamma",lat.count,count3)]][i]<-1/count*cov.table.sem[[paste0("gamma",lat.count,count3)]][i]
      }
  }
}

