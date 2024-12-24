library(Matrix)
library(fda)
library(doParallel)
library(foreach)

numCores<-detectCores()-1
c1<-makeCluster(numCores)
registerDoParallel(c1)


params.estimated.se <-results
n<-400
d.f.se <- readRDS(file.path("realData", paste0("d.f.se_n", n, ".rds")))

model.dat.se<-fsem(eta~~z1,effectType = "fixed_concurrent")%+%
  fsem(eta~~z2,effectType="concurrent")%+%
  fsem(eta~~z3,effectType="concurrent")%+%
  fsem(eta~~z4,effectType="concurrent")%+%
  fsem(eta~~z5,effectType="concurrent")%+%
  fsem(eta~-1+gender,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+cancer,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+diabet,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+heartcon,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+genCancer,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+genDiabet,effectType="linear",scalar.covariate = TRUE)%+%
  fsem(eta~-1+genHeart,effectType="linear",scalar.covariate = TRUE)

dat_list<-list()
dat_list[[1]]<-params.estimated.se

n<-length(dat_list)
dat_list<-list()
dat_list[[1]]<-list()
dat_list[[1]]$model.fit<-model.dat.se
dat_list[[1]]$model.sim<-model.dat.se
dat_list[[1]]$estimation<-params.estimated.se
dat_list[[1]]$simmulation<-d.f.se

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
for (n.out in 1:n) {
  re<-dat_list[[n.out]]
  no.f<-length(re$model.fit$mod$factorModel)
  no.r<-length(re$model.fit$mod$regression)
  lat2<-re$model.fit$var$latents
  lat2.cov<-re$model.fit$var$observed
  n.b<-length(re$estimation$result$params.estimated$intercept[[1]])
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(0,1))
  times<-(seq(0,1,length.out=200))[ind]
  eval1<-t(eval.basis(times,basis))
  samp<-unique(re$simmulation$data$.id)
  n.samp<-length(samp)
  n.t<-length(unique(re$simmulation$data$.t))
  
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
  foreach(k=1:boot, .packages = c("Matrix","fda")) %dopar%{
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
              g[[m]][[i1]]<-re$estimation$hist$hist.sem[[lat.count]]$g.x[[i1]]
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
            count3<-which(lat2.cov==re$model.sim$mod$regression[[j]]$covariate[[jj]])
            d.sem1<-t(eval1)%*%sig.sem[[k]][[lat.count]]%*%eval1  
            if(re$model.sim$mod$regression[[j]]$effect[[jj]]%in%c("concurrent","linear")){
              ev<-eval1
              ga<-gamma2.1[(count+1):(count+n.b)]
              count<-count+n.b
              gamma1[[k]][[lat.count]][[count3]]<-diag(diag(d.sem1)^(-1/2))%*%t(ev)%*%ga
            }else{
              evv_ind<-t(eval.basis((seq(0,times[5],length.out=200))[ind],basis))
              ev<-diag(d.sem1)^(-1/2)*t(eval.basis(times,basis))
              evv<-(diag(d.sem1)^(-1/2))[5]*kronecker(ev[,5],evv_ind)
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
  
  
  stopCluster(c1)
  ress<-readRDS("D:\\Github\\Output_realData\\correctcode\\res_parallel_con_400.rds")
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
        count3<-which(lat2.cov==re$model.sim$mod$regression[[j]]$covariate[[jj]])
        gam.sem[[n.out]][[lat.count]][[count3]]<-list()
        cov.sem[[n.out]][[lat.count]][[count3]]<-list()
        per.sem[[n.out]][[lat.count]][[count3]]<-list()
        res<-lat2[lat.count]
        cov<-lat2.cov[count3]
        gamma1.1<-t(sapply(1:boot,function(k){ress[[k]]$gamma1[[lat.count]][[count3]]}))
        if(re$model.sim$mod$regression[[j]]$effect[[jj]]%in%c("concurrent","linear")){
          gam.sem[[n.out]][[lat.count]][[count3]]<-re$estimation$result$params.estimated$gamma.param.x.std[[lat.count]][[count3]][ind]  
          cov.sem[[n.out]][[lat.count]][[count3]]<-cov(gamma1.1)
        }else{
          gam.sem[[n.out]][[lat.count]][[count3]]<-re$estimation$result$params.estimated$gamma.param.x.std[[lat.count]][[count3]][ind,ind][,5]
          cov.sem[[n.out]][[lat.count]][[count3]]<-cov(gamma1.1)
        }
        gam.sem[[n.out]][[lat.count]][[count3]]<-as.vector(gam.sem[[n.out]][[lat.count]][[count3]])
        cov.sem[[n.out]][[lat.count]][[count3]]<-as.matrix(cov.sem[[n.out]][[lat.count]][[count3]])
        fregion.bands.sem[[n.out]][[lat.count]][[count3]] = fregion::fregion.band(gam.sem[[n.out]][[lat.count]][[count3]],cov.sem[[n.out]][[lat.count]][[count3]], type = c("BEc"),conf.level=c(0.95))
      } 
    }
  }
  for (j in 1:no.f) {
    lambda1.1<-as.matrix(t(sapply(1:boot,function(k){ress[[k]]$lambda1[[j]]})))
    count3<-which(lat2==re$model.fit$mod$factorModel[[j]]$factor)
    res<-re$model.fit$var$indicators[j]
    cov<-lat2[count3]
    if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){
      lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind]
      cov.fac[[n.out]][[j]]<-cov(lambda1.1)
    }else{
      if(re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
        lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind]  
        cov.fac[[n.out]][[j]]<-cov(lambda1.1)
      }else{
        lam.fac[[n.out]][[j]]<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind,ind][,5]
        cov.fac[[n.out]][[j]]<-cov(lambda1.1)
      } 
    }
    lam.fac[[n.out]][[j]]<-as.vector(lam.fac[[n.out]][[j]])
    cov.fac[[n.out]][[j]]<-as.matrix(cov.fac[[n.out]][[j]])
    fregion.bands.fac[[n.out]][[j]] = fregion::fregion.band(lam.fac[[n.out]][[j]],cov.fac[[n.out]][[j]],type = c("BEc"),conf.level=c(0.95))
  }     
}
