###########################################################################
##
## file: LESN.R  (originally this file was called AdHoc.R, but it was
##                formally renamed on Mar 21, 2003)
##
## Copyright (C) 2003     Ben Bolstad
##
## "AdHoc" background correction routines - The LESN routines
##
## Feb 14, 2003 - give reasonable defaults to functions
## Feb 27, 2003 - clean out two commented out functions
## Mar 21, 2003 - Rename file LESN.R to keep consistant naming across methods
##
###########################################################################

bg.correct.shift <- function(object,baseline=0.25){

  objsize <- dim(pm(object))
  
  pm(object) <- matrix(.C("R_shift_down",as.double(as.vector(pm(object))),as.double(baseline),as.integer(objsize[1]),as.integer(objsize[2]))[[1]],objsize[1],objsize[2])
  object
}


bg.correct.stretch <- function(object,baseline=0.25,type=c("linear","exponential","loglinear","logexponential","loggaussian"),theta){
  
  objsize <- dim(pm(object))
  
  if (type == "linear"){
    ty <- 1
  } else if (type == "exponential"){
    ty <- 2
  } else if (type == "loglinear"){
    ty <- 3
  } else if (type=="logexponential"){
    ty <- 4
  } else {
    ty <- 5
    
  }
  
  pm(object) <- matrix(.C("R_stretch_down",as.double(as.vector(pm(object))),as.double(baseline),as.integer(objsize[1]),as.integer(objsize[2]),as.integer(ty),as.double(theta))[[1]],objsize[1],objsize[2])
  object
}




bg.correct.LESN <- function(object,method = 2,baseline = 0.25, theta=4){
  if (method == 2){
    bg.correct.stretch(object, baseline, type="loggaussian",theta=2*theta^2) 
  } else if (method == 1){
    bg.correct.stretch(object, baseline, type="logexponential",theta) 
  } else {
    bg.shift(object,baseline)
  }
}

