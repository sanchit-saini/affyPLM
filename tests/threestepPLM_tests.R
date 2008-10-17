if (.Platform$OS.type != "windows"){
library(affyPLM)

# test threestep and threestepPLM to see if they agree


check.coefs <- function(Pset,Pset2){
  if (any(abs(coefs(Pset) - exprs(Pset2)) > 1e-14)){
   stop("No agreement between threestepPLM and threestep in coefs")
  }	
}

check.resids <- function(Pset,Pset2){
 if (any(resid(Pset) != resid(Pset2))){
   stop("No agreement between threestepPLM and rmaPLM/threestep in residuals")
}

}



data(Dilution)

Pset <- threestepPLM(Dilution)
Pset2 <- threestep(Dilution)
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="tukey.biweight")
Pset2 <- threestep(Dilution,summary.method="tukey.biweight")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="average.log")
Pset2 <- threestep(Dilution,summary.method="average.log")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="rlm")
Pset2 <- threestep(Dilution,summary.method="rlm")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="log.average")
Pset2 <- threestep(Dilution,summary.method="log.average")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="log.median")
Pset2 <- threestep(Dilution,summary.method="log.median")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="median.log")
Pset2 <- threestep(Dilution,summary.method="median.log")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="log.2nd.largest")
Pset2 <- threestep(Dilution,summary.method="log.2nd.largest")
check.coefs(Pset,Pset2)

Pset <- threestepPLM(Dilution,summary.method="lm")
Pset2 <- threestep(Dilution,summary.method="lm")
check.coefs(Pset,Pset2)

#check if threestepPLM agrees with rmaPLM
Pset <- threestepPLM(Dilution)
Pset2 <- rmaPLM(Dilution)

if (any(coefs(Pset) != coefs(Pset2))){
   stop("No agreement between threestepPLM and rmaPLM in coefs")
}


if (any(resid(Pset)[[1]] != resid(Pset2)[[1]])){
   stop("No agreement between threestepPLM and rmaPLM in residuals")
}
}
