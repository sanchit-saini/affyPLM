#test the preprocessing functionality

library(affyPLM)
library(affydata)
data(Dilution)


### NO LONGER SUPPORTED eset <- threestep(Dilution,background.method="RMA.1")
eset <- threestep(Dilution,background.method="RMA.2")
eset <- threestep(Dilution,background.method="IdealMM")
eset <- threestep(Dilution,background.method="MAS")
eset <- threestep(Dilution,background.method="MASIM")
eset <- threestep(Dilution,background.method="LESN2")
eset <- threestep(Dilution,background.method="LESN1")
eset <- threestep(Dilution,background.method="LESN0")

eset <- threestep(Dilution,normalize.method="quantile",background=FALSE)
eset <- threestep(Dilution,normalize.method="quantile.probeset",background=FALSE)
eset <- threestep(Dilution,normalize.method="scaling",background=FALSE)


