#test the preprocessing functionality

library(affyPLM)
data(affybatch.example)


eset <- threestep(affybatch.example,background.method="RMA.1")
eset <- threestep(affybatch.example,background.method="RMA.2")
eset <- threestep(affybatch.example,background.method="IdealMM")
eset <- threestep(affybatch.example,background.method="MAS")
eset <- threestep(affybatch.example,background.method="MASIM")
eset <- threestep(affybatch.example,background.method="LESN2")
eset <- threestep(affybatch.example,background.method="LESN1")
eset <- threestep(affybatch.example,background.method="LESN0")

eset <- threestep(affybatch.example,normalize.method="quantile",background=FALSE)
eset <- threestep(affybatch.example,normalize.method="quantile.probeset",background=FALSE)
eset <- threestep(affybatch.example,normalize.method="scaling",background=FALSE)


