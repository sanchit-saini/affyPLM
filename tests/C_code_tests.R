####
#### This code is messy, possibly incomplete and only for
#### the use of developers.
####
####

test.c.code <-  FALSE
test.PLM.modelmatrix <- FALSE
test.rlm <- FALSE

if (test.c.code){
  
  library(affyPLM)
  narrays <- 10
  nprobes <- 11
  nprobetypes <- 2
  ncols <- 10
  
  MMs <- rnorm(narrays*nprobes*nprobetypes)
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  
                                        #test making intercept column
  matrix(.C("R_PLM_matrix_intercept",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),0)[[1]],ncol=10)
  
                                        #test making an MM covariate column
  matrix(.C("R_PLM_matrix_MM",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.double(MMs))[[1]],ncol=10)
  
                                        # sample effect aka chip effect, aka expression values
  matrix(.C("R_PLM_matrix_sample_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0))[[1]],ncol=10)
  matrix(.C("R_PLM_matrix_sample_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1))[[1]],ncol=10)
  matrix(.C("R_PLM_matrix_sample_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1))[[1]],ncol=10)
  
  
  
                                        #probe-type parameter overall
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(0),integer(narrays),as.integer(1))[[1]],ncol=10)
  
                                        #probe-type parameter within sample
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(1),integer(narrays),as.integer(1))[[1]],ncol=10)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(1),integer(narrays),as.integer(1))[[1]],ncol=10)
  ncols <- 20
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(1),integer(narrays),as.integer(0))[[1]],ncol=20)
  
  
                                        #probe-type-parameter within a chip-level factor (eg treatment, or genotype variable)
  trt.cov <- rep(0:1,5)
  ncols <- 10
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=10)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=10)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=10)
  
  trt.cov <- rep(0:4,2)
  matrix(.C("R_PLM_matrix_probe_type_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(2),as.integer(trt.cov),as.integer(4))[[1]],ncol=10)
  
  
  
  
                                        #probe effects - overall
  ncols <- 11
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  
  
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(0),as.integer(trt.cov),as.integer(4))[[1]],ncol=11)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(0),as.integer(trt.cov),as.integer(4))[[1]],ncol=11)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(0),as.integer(trt.cov),as.integer(4))[[1]],ncol=11)
  
  
  
                                        #probe effects within treatment or genotype factor
  trt.cov <- rep(0:1,5)
  ncols <- 22
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(2),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  
  
                                        #probe effects within probetype
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(3),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(3),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(3),as.integer(trt.cov),as.integer(1))[[1]],ncol=22)
  
  
                                        #probe effects within probetype within treatment or genotype factor variable
  trt.cov <- rep(0:1,5)
  ncols <- 44
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  
  nprobetypes <- 1
  trt.cov <- rep(0:1,5)
  ncols <- 44
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  matrix(.C("R_PLM_matrix_probe_effect",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(-1),as.integer(4),as.integer(trt.cov),as.integer(1))[[1]],ncol=44)
  
  
                                        # copy across chip level variables into model matrix
  nprobetypes <- 1
  trt.cov <- rep(0:1,5)
  ncols <- 10
  X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
  trt.variables <- rnorm(10)
  
  matrix(.C("R_PLM_matrix_chiplevel",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.double(trt.variables),as.integer(1))[[1]],ncol=10)


###
### Build a few design matrices and compare with R model.matrix
###


  for (nprobetypes in 1:2){
    for (narrays in 2:15){
      for (nprobes in 2:20){
        for (constraint.type in c("contr.sum","contr.treatment")){
          if (constraint.type == "contr.sum"){
            ct.type <- -1
          } else {
            ct.type <- 1
          }

          
          ncols <- nprobes -1 + narrays
          X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
          
          X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(1),as.integer(0),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
          
          probe.effect <- factor(rep(1:nprobes,narrays*nprobetypes))
          sample.effect <- factor(rep(rep(c(1:narrays),rep(nprobes,narrays)),nprobetypes))
          if (nprobetypes == 2){
            probe.type.effect <- factor(rep(1:2,c(narrays*nprobes,narrays*nprobes)))
          } else {
            probe.type.effect <- factor(rep(1,narrays*nprobes))
          }
          
          if (any(X!=model.matrix(~ C(sample.effect,constraint.type) + C(probe.effect,constraint.type)))){
            stop("Model matrix function problem ",narrays," ", nprobes, " ", nprobetypes)
          }
          
          X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
          X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(0),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
          
          if (any(X!=model.matrix(~ -1 + C(sample.effect,constraint.type) + C(probe.effect,constraint.type)))){
            stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
          }
          
          ncols <- nprobes
          X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
          X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(0),as.integer(0),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
          
          if (any(X!=model.matrix(~  C(probe.effect,constraint.type)))){
            stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
          }
          
          
          ncols <- nprobes
          X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
          X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(0),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
          
          if (any(X!=model.matrix(~-1+  C(probe.effect,constraint.type)))){
            stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
          }
          
        }
      }
    }
  }
  
###
### Build a few more design matrices and compare with R model.matrix
###

  
  for (narrays in 2:15){
    for (nprobes in 2:20){
      for (constraint.type in c("contr.sum","contr.treatment")){
        probe.effect <- factor(rep(1:nprobes,narrays*nprobetypes))
        sample.effect <- factor(rep(rep(c(1:narrays),rep(nprobes,narrays)),nprobetypes))
        if (constraint.type == "contr.sum"){
          ct.type <- -1
        } else {
          ct.type <- 1
        }

       
        if (nprobetypes == 2){
          probe.type.effect <- factor(rep(1:2,c(narrays*nprobes,narrays*nprobes)))
        } else {
          probe.type.effect <- factor(rep(1,narrays*nprobes))
        }
        ncols <- nprobetypes + nprobes -1
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(1),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~-1+ C(probe.type.effect,constraint.type) +  C(probe.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }
        
        ncols <- nprobetypes + nprobes -1
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(0),as.integer(1),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ C(probe.type.effect,constraint.type) +  C(probe.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }
        
        ncols <- narrays + nprobetypes + nprobes -2
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(1),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ -1 + C(sample.effect,constraint.type) + C(probe.type.effect,constraint.type) +  C(probe.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }
        
        ncols <- narrays + nprobetypes + nprobes -2
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(1),as.integer(1),as.integer(1),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ + C(sample.effect,constraint.type) + C(probe.type.effect,constraint.type) +  C(probe.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }
        
        ncols <- narrays + nprobetypes -1
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(1),as.integer(1),as.integer(0),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ C(sample.effect,constraint.type) + C(probe.type.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }    
        
        ncols <- narrays + nprobetypes -1
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(1),as.integer(1),as.integer(0),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ -1 + C(sample.effect,constraint.type) + C(probe.type.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }    
        
        ncols <- nprobetypes
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(1),as.integer(0),as.integer(1),as.integer(0),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~ C(probe.type.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }    
      
        ncols <- nprobetypes
        X <- matrix(0,narrays*nprobes*nprobetypes,ncols)
        X <- matrix(.C("R_PLM_Matrix_constructtest",as.double(as.vector(X)), as.integer(narrays),as.integer(nprobes),as.integer(nprobetypes),as.integer(0),as.integer(0),as.integer(1),as.integer(0),as.integer(ct.type))[[1]],ncol=ncols)
        
        if (any(X!=model.matrix(~-1 +  C(probe.type.effect,constraint.type)))){
          stop("Model matrix function problem ",narrays," ", nprobes," ", nprobetypes)
        }    
        
      }
      
      
    }
  }


  narrays <- 2
  nprobes <- 7
  nprobetypes <- 2
  
  probe.effect <- factor(rep(1:nprobes,narrays*nprobetypes))
  sample.effect <- factor(rep(rep(c(1:narrays),rep(nprobes,narrays)),nprobetypes))
  if (constraint.type == "contr.sum"){
    ct.type <- -1
  } else {
    ct.type <- 1
  }
  
  
  if (nprobetypes == 2){
    probe.type.effect <- factor(rep(1:2,c(narrays*nprobes,narrays*nprobes)))
  } else {
    probe.type.effect <- factor(rep(1,narrays*nprobes))
  }
  

  model.matrix(~-1 +probe.effect/probe.type.effect)


  library(affyPLM)
  output <- verify.output.param(list(weights = FALSE, residuals = FALSE, varcov = "none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  data(Dilution)

  # fit a PM ~ samples model
  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0), probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)


  sample.effect <- rep(1:4,c(16,16,16,16))
  probe.effect <- rep(1:16,4)

  library(MASS)
  fit <- rlm(as.vector(log2(pm(Dilution)[1:16,])) ~ -1 + factor(sample.effect))
  
  if (any(Fitresults[[1]][1,] != coef(fit))){
    stop("Problem in model fitting procedure")
  }

  sample.effect <- rep(1:4,c(20,20,20,20))
  probe.effect <- rep(1:20,4)
  fit <- rlm(as.vector(log2(pm(Dilution)[201781:201800,])) ~ -1 + factor(sample.effect))
   if (any(Fitresults[[1]][12625,] != coef(fit))){
    stop("Problem in model fitting procedure")
  }


  # fit a samples + probes model
  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0), probe.type.levels=list(),probe.trt.levels=list())

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  sample.effect <- rep(1:4,c(16,16,16,16))
  probe.effect <- rep(1:16,4)

  fit <- rlm(as.vector(log2(pm(Dilution)[1:16,])) ~ -1 + factor(sample.effect)+C(factor(probe.effect),"contr.sum"))
  
  if (any(abs(Fitresults[[1]][1,] -coef(fit)[1:4]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

   if (any(abs(as.vector(Fitresults[[2]][[1]]) - coef(fit)[5:19]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }



  
  sample.effect <- rep(1:4,c(20,20,20,20))
  probe.effect <- rep(1:20,4)
  fit <- rlm(as.vector(log2(pm(Dilution)[201781:201800,])) ~ -1 + factor(sample.effect)+C(factor(probe.effect),"contr.sum"))
  
  if (any(abs(Fitresults[[1]][12625,]  -coef(fit)[1:4])> 1e-13)){
    stop("Problem in model fitting procedure")
  }

  if (any(abs(as.vector(Fitresults[[2]][[12625]])- coef(fit)[5:23])>1e-13)){
    stop("Problem in model fitting procedure")
  }

  # fit an MM ~ samples model
  R.model <- list(mmorpm.covariate=0,response.variable=-1,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0), probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)


  sample.effect <- rep(1:4,c(16,16,16,16))
  probe.effect <- rep(1:16,4)

  library(MASS)
  fit <- rlm(as.vector(log2(mm(Dilution)[1:16,])) ~ -1 + factor(sample.effect))
  
  if (any(abs(Fitresults[[1]][1,] - coef(fit)) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  sample.effect <- rep(1:4,c(20,20,20,20))
  probe.effect <- rep(1:20,4)
  fit <- rlm(as.vector(log2(mm(Dilution)[201781:201800,])) ~ -1 + factor(sample.effect))
   if (any(abs(Fitresults[[1]][12625,] - coef(fit))>1e-13)){
    stop("Problem in model fitting procedure")
  }

  # fit a MM ~ samples + probes model
  R.model <- list(mmorpm.covariate=0,response.variable=-1,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0), probe.type.levels=list(),probe.trt.levels=list())

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  sample.effect <- rep(1:4,c(16,16,16,16))
  probe.effect <- rep(1:16,4)

  fit <- rlm(as.vector(log2(mm(Dilution)[1:16,])) ~ -1 + factor(sample.effect)+C(factor(probe.effect),"contr.sum"))
  
  if (any(abs(Fitresults[[1]][1,]-coef(fit)[1:4]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

   if (any(abs(as.vector(Fitresults[[2]][[1]])- coef(fit)[5:19])>1e-13)){
    stop("Problem in model fitting procedure")
  }



  
  sample.effect <- rep(1:4,c(20,20,20,20))
  probe.effect <- rep(1:20,4)
  fit <- rlm(as.vector(log2(mm(Dilution)[201781:201800,])) ~ -1 + factor(sample.effect)+C(factor(probe.effect),"contr.sum"))
  
  if (any(abs(Fitresults[[1]][12625,]- coef(fit)[1:4])>1e-13)){
    stop("Problem in model fitting procedure")
  }

  if (any(abs(as.vector(Fitresults[[2]][[12625]])-coef(fit)[5:23])>1e14)){
    stop("Problem in model fitting procedure")
  }


  # a treatment model
  treatment.effect <- c(1,1,2,2)

  covariates <- model.matrix(~ -1 + as.factor(treatment.effect))

    R.model <- list(mmorpm.covariate=0,response.variable=-1,which.parameter.types=as.integer(c(0,1,0,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =covariates, probe.type.levels=list(),probe.trt.levels=list())
  
Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  treatment.effect <- rep(c(1,1,2,2),c(16,16,16,16))
  fit <- rlm(as.vector(log2(mm(Dilution)[1:16,])) ~ -1 + factor(treatment.effect))
  
  if (any(abs(Fitresults[[1]][1,]-coef(fit)[1:2]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

   output <-  verify.output.param(list(weights = FALSE, residuals = FALSE, varcov = "none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)


  # a treatment + probes model with contr.treatment constraint
  R.model <- list(mmorpm.covariate=0,response.variable=-1,which.parameter.types=as.integer(c(0,1,0,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =covariates, probe.type.levels=list(),probe.trt.levels=list())
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  
  treatment.effect <- rep(c(1,1,2,2),c(20,20,20,20))
  probe.effect <- rep(1:20,4)
  fit <- rlm(as.vector(log2(mm(Dilution)[201761:201780,])) ~ -1 + factor(treatment.effect)+C(factor(probe.effect),"contr.treatment"))
  
  if (any(abs(Fitresults[[1]][12624,]-coef(fit)[1:2]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  if (any(abs(as.vector(Fitresults[[2]][[12624]])-coef(fit)[3:21])>1e14)){
    stop("Problem in model fitting procedure")
  }


  
  
  # MM + samples + probes
  R.model <- list(mmorpm.covariate=1,response.variable=1,which.parameter.types=as.integer(c(0,0,0,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0), probe.type.levels=list(),probe.trt.levels=list())


  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  
 sample.effect <- rep(1:4,c(16,16,16,16))
  probe.effect <- rep(1:16,4)

  fit <- rlm(as.vector(log2(pm(Dilution)[1:16,])) ~ -1 + as.vector(log2(mm(Dilution)[1:16,])))
  

  if (any(abs(as.vector(Fitresults[[6]][[1]]) - coef(fit)) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  R.model <- list(mmorpm.covariate=1,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())


  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  fit <- rlm(as.vector(log2(pm(Dilution)[1:16,])) ~ -1 + as.vector(log2(mm(Dilution)[1:16,]))+ as.factor(sample.effect))
  if (any(abs(as.vector(Fitresults[[6]][[1]]) - coef(fit)[1]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

 if (any(abs(as.vector(Fitresults[[1]][1,]) - coef(fit)[2:5]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  R.model <- list(mmorpm.covariate=1,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())


  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  fit <- rlm(as.vector(log2(pm(Dilution)[1:16,])) ~ -1 + as.vector(log2(mm(Dilution)[1:16,]))+ as.factor(sample.effect) + C(as.factor(probe.effect),"contr.sum"))
  if (any(abs(as.vector(Fitresults[[6]][[1]]) - coef(fit)[1]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  if (any(abs(as.vector(Fitresults[[1]][1,]) - coef(fit)[2:5]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }


  ## PM and MM are response

  
  sample.effect <- rep(1:4,c(32,32,32,32))
  probe.effect <- rep(1:16,8)
  

  
  # PMMM ~ -1 + samples
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  
  fit <- rlm(as.vector(rbind(log2(pm(Dilution)[1:16,]),log2(mm(Dilution)[1:16,]))) ~ -1 + as.factor(sample.effect))
  if (any(abs(as.vector(Fitresults[[1]][1,]) - coef(fit)) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  # PMMM ~ -1 + samples +PROBES
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  fit <- rlm(as.vector(rbind(log2(pm(Dilution)[1:16,]),log2(mm(Dilution)[1:16,]))) ~ -1 + as.factor(sample.effect)+C(as.factor(probe.effect),"contr.sum") )
  if (any(abs(as.vector(Fitresults[[1]][1,]) - coef(fit)[1:4]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

  # a probe.type effect
  probe.type.effect <- rep(rep(1:2,c(16,16)),4)

  # PMMM ~ -1 + samples + probe.type + PROBES
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,1,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,-1,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  fit <- rlm(as.vector(rbind(log2(pm(Dilution)[1:16,]),log2(mm(Dilution)[1:16,]))) ~ -1 + as.factor(sample.effect)+ C(as.factor(probe.type.effect),"contr.sum")+ C(as.factor(probe.effect),"contr.sum") )
  
   if (any(abs(as.vector(Fitresults[[1]][1,]) - coef(fit)[1:4]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }

   if (any(abs(as.vector(Fitresults[[6]][1,]) - coef(fit)[5]) > 1e-13)){
    stop("Problem in model fitting procedure")
  }



  #### store weights PM ~ -1 + samples

  output <-  verify.output.param(list(weights = TRUE, residuals = TRUE, varcov = "none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  
  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
 

#### store weights PMMM ~ -1 + samples

  output <-  verify.output.param(list(weights = TRUE, residuals = TRUE, varcov = "none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
#### store weights PMMM ~ -1 + samples + probe.type + probes

  output <- verify.output.param(list(weights = TRUE, residuals = TRUE, varcov ="none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,1,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,1,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  ## PM ~ -1 + treatment + probes in treatment
  output <- verify.output.param(list(weights = TRUE, residuals = TRUE, varcov = "none", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  treatment.effect <- c(1,1,2,2)

  covariates <- model.matrix(~ -1 + as.factor(treatment.effect))

  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,1,0,0,1)),strata=as.integer(c(0,0,0,0,2)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =covariates,probe.type.levels=list(),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  ## PMMM ~ -1 + treatment + probes in treatment
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  treatment.effect <- c(1,1,2,2)

  covariates <- model.matrix(~ -1 + as.factor(treatment.effect))

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,1,0,0,1)),strata=as.integer(c(0,0,0,0,2)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =covariates,probe.type.levels=list(),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
   ## PMMM ~ -1 + treatment + probe.effect in treatment
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  treatment.effect <- c(1,1,2,2)

  covariates <- model.matrix(~ -1 + as.factor(treatment.effect))

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,1,0,1,0)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,-1,0)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =covariates,probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
   ## PMMM ~ -1+ probes
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
   ## PMMM ~ -1+ probes.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,0)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  ## PMMM ~ -1+ probes.type + probes
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=NULL,max.probe.type.trt.factor=0,probe.trt.factor=NULL,max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list(),probe.trt.levels=list())

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

    ## PMMM ~ -1+ probes.type + probes     with both within treatment factor
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,1)),strata=as.integer(c(0,0,0,2,2)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)


    ## PMMM ~ -1+ probes.type + probes     with both within treatment factor and probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,1)),strata=as.integer(c(0,0,0,2,4)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)


  ## PMMM ~ -1+ probes.type + probes     probe.types within treatment factor and probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,1)),strata=as.integer(c(0,0,0,2,3)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

 ## PMMM ~ -1+ probes.type + probes     probe.types within samples and probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,1,1)),strata=as.integer(c(0,0,0,1,3)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)



  ## PMMM ~ intercept + probes.type + probes     probe.types within samples and probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,0,0,1,1)),strata=as.integer(c(0,0,0,1,3)),constraints=as.integer(c(0,0,0,-1,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

    ## PMMM ~ intercept + probes      probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,0,0,0,1)),strata=as.integer(c(0,0,0,0,3)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

     ## PMMM ~ -1+ probes      probes also within probe.type
  output <- list(weights = TRUE, residuals = TRUE, varcov = c("none","chiplevel", "all"), resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,0,0,1)),strata=as.integer(c(0,0,0,0,3)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  # now play with varcov output
  output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

 # now play with varcov output and treatment
  output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
 treatment.effect <- c(1,1,2,2)

  covariates <- model.matrix(~ -1 + as.factor(treatment.effect))

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,1,0,0,0)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =covariates,probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  

  # now play with varcov output and an intercept
  output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,0,1,0,1)),strata=as.integer(c(0,0,0,0,0)),constraints=as.integer(c(0,0,-1,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=1,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))
  
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)



 # now play with varcov output and treatment and intercept
  output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
 treatment.effect <- c(1,1,2,2)

  covariates <- matrix(model.matrix(~ as.factor(treatment.effect))[,2])
  colnames(covariates) <- "trt_2"

  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,1,0,0,0)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =covariates,probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  

 # now play with varcov all option output and treatment and intercept
  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
 treatment.effect <- c(1,1,2,2)
   
  covariates <- matrix(model.matrix(~ as.factor(treatment.effect))[,2])
  colnames(covariates) <- "trt_2"
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,1,0,0,0)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,0)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =covariates,probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

 # now play with varcov all option output and samples and intercept
  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(1,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,-1,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

       # now play with varcov all option output and samples and intercept, MM covarite
    output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=1,response.variable=1,which.parameter.types=as.integer(c(1,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,-1,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)







  ## now play with varcov all option output and samples and intercept, MM covariate and input chip weights
  output <- verify.output.param(list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE))
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=c(1,1,0.5,0.5),weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  ## now play with varcov all option output and samples and intercept, MM covariate and input chip weights
  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=runif(201800))
  
  R.model <- list(mmorpm.covariate=0,response.variable=1,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)



  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=c(rep(c(1,0.5),c(201800,201800))))
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  

  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="cuberoot", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  

  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log10", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  
  R.model <- list(mmorpm.covariate=0,response.variable=0,which.parameter.types=as.integer(c(0,0,1,0,1)),strata=as.integer(c(0,0,0,2,0)),constraints=as.integer(c(0,0,0,0,-1)),probe.type.trt.factor=as.integer(c(0,0,1,1)),max.probe.type.trt.factor=1,probe.trt.factor=as.integer(c(0,0,1,1)),max.probe.trt.factor=0,chipcovariates =matrix(0,0,0),probe.type.levels=list("blah"=c("A","B")),probe.trt.levels=list("blah"=c("A","B")))

  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
}



if (test.PLM.modelmatrix){

  library(affyPLM);data(Dilution)
  
  #PLM.designmatrix3(Dilution)
  
  #PLM.designmatrix3(Dilution,model=MM ~ PM -1 + samples +probe.type:probes)
  
  #PLM.designmatrix3(Dilution,model=MM ~ PM -1 + samples:probe.type + liver:probe.type:probes + liver:samples)
  #PLM.designmatrix3(Dilution,model=MM ~ PM + samples:probe.type + liver:probe.type:probes + liver + samples)




  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  #blah <- c(1,5,5,1)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ probes + blah,constraint.type=c(probes="contr.sum"))
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 + blah:probe.type)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 +probes:probe.type)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 +probes:blah)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 +probes:probe.type:blah)
  #output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
#  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ samples,constraint.type=c(samples="contr.sum"))
#  R.model <-  PLM.designmatrix3(Dilution,model=PMMM ~ blah)
#   R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 + samples)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ probes + blah)
  #R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ probes + blah)
  #Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  library(affyPLM);data(Dilution)
   output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ probe.type + probe.type:probes + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  library(affyPLM);data(Dilution)
  output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ samples:probe.type + probe.type:probes + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  
 output <- list(weights = TRUE, residuals = TRUE, varcov ="chiplevel", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  blah <- c(1,2,2)
  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ blah:probe.type + probe.type:probes + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  
  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  blah <- c(1,2,2)
  R.model <- PLM.designmatrix3(Dilution,model=PM ~ -1 + probes + MM + blah)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)

  
  output <- list(weights = TRUE, residuals = TRUE, varcov ="all", resid.SE = TRUE)
  modelparam <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)
  blah <- c(1,2,2)
  R.model <- PLM.designmatrix3(Dilution,model=PM ~ -1 + probes + MM + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  




#test some of the verification functions


  output <- verify.output.param()
  modelparam <- verify.model.param(Dilution,PM ~ -1 + probes + MM + samples)
  R.model <- PLM.designmatrix3(Dilution,model=PM ~ -1 + probes + MM + samples)
 
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  ##verify.model.param(Dilution,PM ~ -1 + probes + MM + samples,model.param=list(weights.probe=rep(1,10)))

   modelparam <- verify.model.param(Dilution,PMMM ~ -1 + probes + samples,model.param=list(weights.chip=c(1,2,3),weights.probe=rep(1,2400*2)))
  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 + probes + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)


  
 modelparam <- verify.model.param(Dilution,PM ~ -1 + probes + samples,model.param=list())
  R.model <- PLM.designmatrix3(Dilution,model=PM ~ -1 + probes + samples)
  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  ## probes <- rep(1:16,3)
  ##  chips <- rep(1:3,c(16,16,16))

  ## library(MASS)
  
  ##fit <- rlm(log2(as.vector(pm(Dilution,"HG2188-HT2258_at"))) ~ -1 + as.factor(chips) + C(as.factor(probes),"contr.sum"))
  

#test creating a PLMset based on the output from rlm_PLMset

###  x <- new("PLMset")
###  x@chip.coefs=Fitresults[[1]]
###  x@probe.coefs= Fitresults[[2]]
###  x@weights=Fitresults[[3]]
###  x@se.chip.coefs=Fitresults[[4]]
###  x@se.probe.coefs=Fitresults[[5]]
###  x@exprs=Fitresults[[6]]
###  x@se.exprs=Fitresults[[7]]
###  x@residuals=Fitresults[[8]]
###  x@residualSE=Fitresults[[9]]
###  x@varcov = Fitresults[[10]]
###  x@cdfName = Dilution@cdfName
 ### x@phenoData = Dilution@phenoData
 ### x@annotation = Dilution@annotation
###  x@description = Dilution@description
###  x@notes = Dilution@notes
###  x@nrow= Dilution@nrow
###  x@ncol= Dilution@ncol
### x@model.description = c(x@model.description, list(R.model=R.model))
###  image(x)




###  data(Dilution)
###  output <- verify.output.param()
###  modelparam <- verify.model.param(Dilution,PMMM ~ -1 + probe.type:probes + samples + samples:probe.type,model.param=list())
###  R.model <- PLM.designmatrix3(Dilution,model=PMMM ~ -1 + probe.type:probes + samples+ samples:probe.type)
###  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
###  output <- verify.output.param()
###  modelparam <- verify.model.param(Dilution,MM ~ -1 + probes + samples,model.param=list())
###  R.model <- PLM.designmatrix3(Dilution,model=MM ~ -1 + probes + samples)
###  Fitresults <- .Call("rlm_PLMset",pm(Dilution),mm(Dilution),probeNames(Dilution),length(geneNames(Dilution)),R.model,output,modelparam)
  
  

###  x <- new("PLMset")
 ### x@chip.coefs=Fitresults[[1]]
###  x@probe.coefs= Fitresults[[2]]
###  x@weights=Fitresults[[3]]
###  x@se.chip.coefs=Fitresults[[4]]
###  x@se.probe.coefs=Fitresults[[5]]
###  x@exprs=Fitresults[[6]]
###  x@se.exprs=Fitresults[[7]]
 ### x@residuals=Fitresults[[8]]
###  x@residualSE=Fitresults[[9]]
###  x@varcov = Fitresults[[10]]
###  x@cdfName = Dilution@cdfName
###  x@phenoData = Dilution@phenoData
###  x@annotation = Dilution@annotation
###  x@description = Dilution@description
###  x@notes = Dilution@notes
###  x@nrow= Dilution@nrow
###  x@ncol= Dilution@ncol
###  x@model.description = c(x@model.description, list(R.model=R.model))
###  image(x)
###  image(x,type="pos.resids")
###  image(x,type="neg.resids")
###  image(x,type="sign.resids")

###  resid(x,"1091_at")



###  weights(x,c("1091_at","1092_at"))


###  image(x,type="resids",standardize=TRUE)




  
  
}





if (test.rlm){


  library(affyPLM);data(Dilution)

  y <- as.vector(log2(pm(Dilution)[1:16,]))

  w <- runif(64)

  probes <- rep(1:16,4)
  samples <- rep(1:4,c(16,16,16,16))

  x <- model.matrix( ~ -1 + as.factor(samples) + C(as.factor(probes),"contr.sum"))
  x <- as.vector(x)

  cols <- 19
  rows <- 64
  

#  rlm_wfit_R(double *x, double *y, double *w, int *rows, int *cols, double *out_beta, double *out_resids, double *out_weights)

  fit1 <-.C("rlm_wfit_R",as.double(x),as.double(y),as.double(w),as.integer(rows),as.integer(cols),double(cols),double(rows),double(rows))


  library(MASS)

  fit2 <- rlm(y ~ -1 + as.factor(samples) + C(as.factor(probes),"contr.sum"),weights=w,wt.method="case")

  if (any(abs(coef(fit2) - fit1[[6]]) > 10e-14)){
    stop("Weighted RLM did not work")
  }





  y <- as.vector(log2(pm(Dilution,"1001_at")))
  x <- as.vector(log2(mm(Dilution,"1001_at")))

  rlm(y ~ -1 + x + as.factor(samples) + C(as.factor(probes),"contr.sum"))










  
}



