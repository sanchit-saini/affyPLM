bg.correct.gcrma<-  function(object,affinity.info=NULL,
                            affinity.source=c("reference","local"),
                            NCprobe=NULL,
                            type=c("fullmodel","affinities","mm","constant"),
                            k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
                            GSB.adjust=TRUE,rho=0.7,
                            optical.correct=TRUE,verbose=TRUE,fast=TRUE)
{

  bg.adjust.gcrma(object,
                            affinity.info=affinity.info,
                            affinity.source=affinity.source, NCprobe=NCprobe,
                            type=type,k=k,stretch=stretch,correction=correction,
                            GSB.adjust=GSB.adjust,rho=rho,
                            optical.correct=optical.correct,verbose=verbose,fast=fast)


  
  
}
