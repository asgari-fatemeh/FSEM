library(Matrix)
library(lavaan)
library(dplyr)

output<-results
d.f.se2=readRDS("realData\\SEff_matching_n400.rds")

years<-seq(2008,2022,by=2)
chi<-rep(NA,length(years))
rmsea<-chi
srmr<-chi
cfi<-chi
ifi<-chi
rni<-cfi
gfi<-cfi
tli<-cfi
for (ii in 1:length(years)) {
  t=years[ii]
  d<-sapply(1:5, function(j){d.f.se2[[paste0("p",j,"_",t)]]})
  x<-list()
  x$x1<-d.f.se2$gender
  x$x2<-d.f.se2$c_2008
  x$x3<-d.f.se2$d_2008
  x$x4<-d.f.se2$h_2008
  x$x5<-x$x1*x$x2
  x$x6<-x$x1*x$x3
  x$x7<-x$x1*x$x4
  xx<-sapply(1:7, function(i){x[[i]]})
  #observed covariance matrix
  ob<-cbind(d,xx)
  #S<-cov(ob,use="pairwise.complete.obs")
  S<-cov(na.omit(ob))
  
  #predicted covariance matrix
  cov_eps<-diag(sapply(1:5,function(j){output$result$params.estimated$totparam$sigma.error[[j]]}))
  
  times<-seq(2008, 2022, length.out = 200)
  cov_vareps<-diag(sapply(1:5,function(i){lapply(1:5,function(j){approx(times, diag(output$result$params.estimated.eval$ker.fac.std[[paste0("z",j)]]), xout = years)$y})[[i]][ii]}))
  lam<-sapply(1:5,function(i){lapply(1:5,function(j){approx(times,output$result$params.estimated.eval$coef.fac.std[[paste0("z",j)]]$eta, xout = years)$y})[[i]][ii]})
  gam<-sapply(1:7,function(i){lapply(1:7,function(j){approx(times,output$result$params.estimated.eval$coef.sem.std$eta[[j]], xout = years)$y})[[i]][ii]})
  
  
  
  cov_zeta<-approx(times, diag(output$result$params.estimated.eval$ker.sem.std$eta), xout = years)$y[ii]
  cov.x<-diag(S[(6:12),(6:12)])
  cov_et<-sum(gam^2*cov.x)+cov_zeta
  cov_yx<-lam%*%t(gam*cov.x)
  cov_xx<-S[(6:12),(6:12)]
  cov_yy1<-(lam%*%t(lam))*cov_et
  cov_yy<-cov_yy1+cov_vareps+cov_eps
  c1<-cbind(cov_yy,cov_yx)
  c2<-cbind(t(cov_yx),cov_xx)
  
  Sigma<-rbind(c1,c2)
  
  p<-ncol(S)
  k<-45   #number of estimated parameters
  df<-p*(p+1)/2-k
  
  #chi-square
  F_ML<-log(det(Sigma))+sum(diag(S %*% solve(Sigma)))-log(det(S))-ncol(S)
  N<-dim(na.omit(ob))[1]
  chi_square<-((N-1)*F_ML)
  chi[ii]<-chi_square/df
  #less than 3 is good
  
  #Root Mean Square Error of Approximation (RMSEA)
  rmsea[ii]<-sqrt((chi_square-df)/(df*(N-1)))
  #less than 0.05 is good
  
  #Standardized Root Mean Square Residual (SRMR)
  srmr[ii]<-sqrt(sum((S-Sigma)^2)/(p*(p+1)/2))
  #less than 0.08 is good
  
  gfi[[ii]]<-1-sum((S-Sigma)^2)/sum(S^2)
    
  #null model
  SS<-S
  SSigma<-diag(diag(S))
  #SSigma[(6:12),(6:12)]<-S[(6:12),(6:12)]
  FF_ML<-log(det(SSigma))+sum(diag(SS %*% solve(SSigma)))-log(det(SS))-ncol(SS)
  chi_square_null<-((N-1)*FF_ML)
  p<-ncol(SS)
  k<-5+7*8/2
  df_null<-p*(p+1)/2-k

  cfi[[ii]]<-1-(chi_square-df)/(chi_square_null-df_null)
  ifi[[ii]]<-(chi_square_null-chi_square)/(chi_square_null-df_null)
  tli[[ii]]<-((chi_square_null/df_null)-(chi_square/df))/((chi_square_null/df_null)-1)
}
chi_m<-sum(chi)/length(years)
rmsea_m<-sum(rmsea)/length(years)
srmr_m<-sum(srmr)/length(years)
cfi_m<-sum(cfi)/length(years)
ifi_m<-sum(ifi)/length(years)
gfi_m<-sum(gfi)/length(years)
tli_m<-sum(tli)/length(years)
print(round(c(chi_m,rmsea_m,srmr_m,cfi_m,ifi_m,gfi_m,tli_m),digits=3))  

plot(years,ifi,"l",ylim=c(0.7,1),xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
abline(h=0.9, col="red", lty=2,lwd=3)

plot(years,srmr,"l",ylim=c(0,0.09),xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
abline(h=0.08, col="red", lty=2,lwd=3)

###lavaan fitting model
data<-d.f.se2%>%dplyr::select(p1_2008,p2_2008,p3_2008,p4_2008,p5_2008,gender,c_2008,d_2008,h_2008)
data$gc<-data$gender*data$c_2008
data$gd<-data$gender*data$d_2008
data$gh<-data$gender*data$h_2008

model<-'
et=~p1_2008+p2_2008+p3_2008+p4_2008+p5_2008
et~gender+c_2008+d_2008+h_2008+gc+gd+gh
'
fit<-sem(model,data=data)
inds<-fitMeasures(fit,c("chisq","df","pvalue","cfi","ifi","rni","rmsea","srmr","gfi","tli"))


s.hat = fit@Fit@Sigma.hat[[1]]
s.obs = fit@SampleStats@cov[[1]] #cov(data[complete.cases(data),])
s.null = rbind( cbind(diag(diag(s.obs[1:5,1:5])),matrix(0,5,7)),
                cbind(matrix(0,7,5), s.obs[6:12,6:12])) 


F = log(det(s.hat))-log(det(s.obs))+sum(diag(s.obs%*%solve(s.hat)))-ncol(s.hat)
F.null = log(det(s.null))-log(det(s.obs))+sum(diag(s.obs%*%solve(s.null)))-ncol(s.null)
chi2 = (fit@SampleStats@nobs[[1]] ) * F
chi2.null = (fit@SampleStats@nobs[[1]] ) * F.null

gfi = 1-F/F.null

gfi.lav=1-sum((s.obs-s.hat)^2)/sum(s.obs^2)

F.gls = 0.5*sum(diag(((s.obs-s.hat)%*%solve(s.obs))%*%((s.obs-s.hat)%*%solve(s.obs))))
F.gls.null = 0.5*sum(diag(((s.obs-s.null)%*%solve(s.obs))%*%((s.obs-s.null)%*%solve(s.obs))))


gfi.lav=1-F.gls/F.gls.null


#null
s.hat.null=rbind(cbind(diag(diag(s.obs[(1:5),(1:5)])),matrix(rep(0,35),5,7)),cbind(matrix(rep(0,35),7,5),s.obs[(6:12),(6:12)]))
s.obs.null=s.obs
F.null = log(det(s.hat.null))-log(det(s.obs.null))+sum(diag(s.obs.null%*%solve(s.hat.null)))-ncol(s.hat.null)
chi2.null = (fit@SampleStats@nobs[[1]] ) * F.null
cfi.lav<-1-(chi2-df)/(chi2.null-df_null)
ifi.lav<-(chi2.null-chi2)/(chi2.null-df_null)
rni.lav<-1-(chi2-df)/(chi2.null-df_null)
