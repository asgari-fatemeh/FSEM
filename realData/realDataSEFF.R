library(dplyr)
library(fda)
library(Matrix)
library(mvtnorm)
library(tmvtnorm)
library(MASS)

dir.path = ifelse(grepl("^UIO",Sys.info()["nodename"], perl = TRUE), "realData")

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")


n<-400
d.f.se = readRDS(file.path(dir.path, paste0("d.f.se_n",n,".rds")))

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
  


initial.parameter.se<-initial.param(model = model.dat.se,data=d.f.se$data,n.b=6)
params.estimated.se<-em.estimation(model=model.dat.se,data=d.f.se$data,x.data=d.f.se$covariate,
                                   initial.parameter=initial.parameter.se,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular.missing")
results<-params.estimated.se
