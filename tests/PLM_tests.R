# this file tests fitPLM and the PLMset object

library(affyPLM)

data(affybatch.example)


Pset <- fitPLM(affybatch.example)

#check accessors for parameters and se

coefs(Pset)[1:5,]
se(Pset)[1:5,]
coefs.probe(Pset)[1:5]
se.probe(Pset)[1:5]
coefs.const(Pset)
se.const(Pset)

#accessors for weights and residuals

weights(Pset)[[1]][1:5,]
resid(Pset)[[1]][1:5,]


#test varcov

Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,output.param=list(varcov="chiplevel"))
varcov(Pset)[1:3]


#test each of the possible weight functions
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Huber"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="fair"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Cauchy"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Geman-McClure"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Welsch"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Tukey"))
Pset <- fitPLM(affybatch.example,background=FALSE,normalize=FALSE,model.param=list(psi.type="Andrews"))

# a larger example to do some testing of the graphical functions

data(Dilution)

Pset <- fitPLM(Dilution)

#testing the image capabilities

image(Pset,which=2)
image(Pset,which=2,type="resids")
image(Pset,which=2,type="pos.resids")
image(Pset,which=2,type="neg.resids")
image(Pset,which=2,type="resids",use.log=FALSE,add.legend=TRUE)

boxplot(Pset)
Mbox(Pset)


#test some non-default models functions
# no preprocessing for speed

Pset <- fitPLM(Dilution, PM ~ -1 + probes + liver,background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]


Pset <- fitPLM(Dilution, PM ~ -1 + probes + liver + scanner,background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]

#checking the constraints
Pset <- fitPLM(Dilution, PM ~ -1 + probes + liver + scanner,constraint.type=c(default="contr.sum"),background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]

Pset <- fitPLM(Dilution, PM ~ -1 + liver + scanner,constraint.type=c(default="contr.sum"),background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]
coefs.probe(Pset) # should be empty

Pset <- fitPLM(Dilution, PM ~ -1 + liver + scanner,constraint.type=c(probes="contr.treatment"),background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]
coefs.probe(Pset) # should be empty


Pset <- fitPLM(Dilution, PM ~ -1 + probes + liver + scanner,constraint.type=c(probes="contr.treatment"),background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]
coefs.probe(Pset)[1:16] 


scanner2 <- c(1,2,1,2) 
Pset <- fitPLM(Dilution, PM ~ -1 + probes + liver + scanner2,constraint.type=c(probes="contr.sum"),background=FALSE,normalize=FALSE)
coefs(Pset)[1:5,]
se(Pset)[1:5,]
coefs.probe(Pset)[1:16] 

#
#Pset <- fitPLM(Dilution,model=PM~-1+probes+scanner,normalize=FALSE,background=FALSE,model.param=list(se.type=3))
#se(Pset)[1:10,]

#check that fitPLM rlm agrees with threestep rlm and threestepPLM rlm


Pset <- fitPLM(affybatch.example)
eset <- threestep(affybatch.example,summary.method="rlm")
Pset2 <- threestepPLM(affybatch.example,summary.method="rlm")

if (any(abs(coefs(Pset) - exprs(eset)) > 1e-14)){
  stop("no agreement between fitPLM and threestep")
}

if (any(abs(coefs(Pset) - coefs(Pset2)) > 1e-14)){
  stop("no agreement between fitPLM and threestep")
}
