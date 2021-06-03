library("QuantumOps")
#problem,graph,Error,error,e_mode,tidx  <depth>, <RC>, <extraCycles>, <state>, <samples>

state <- ""
depth <- ""
RC <- FALSE
extraCycles <- FALSE
samples <- 1

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 6)
	print("Got less than 6 arguments")
problem <- args[1]
graph <- args[2]
Error <-args[3]
error <- as.numeric(args[4])
e_mode <- as.integer(args[5])
tidx <- as.integer(args[6])
if(length(args) >= 7)
	depth <- as.integer(args[7])
if(length(args) >= 8)
	RC <- as.integer(args[8]) > 0
if(length(args) >= 9)
	extraCycles <- as.integer(args[9]) > 0
if(length(args) >= 10)
	state <- args[10]
if(length(args) >= 11)
	samples <- as.integer(args[11])

#Coherent noise goes from 0 - error*pi
if(Error == 'coherent' | Error == 'coherentx')
	error <- error * pi

#100 data points between 0 and error
nPoints <- 100			
if(e_mode == 4)
	error <- seq(0,error,length.out=nPoints)[tidx]

#Problem specific setup
Exp <- FALSE
#If problem is anything other than QAOA, graph serves role as nQubits
if(problem != 'qaoa' && problem != 'qaoaDecomposed'){
    graph <- as.integer(graph)
}else{  #Also find Expectation Value as metric when running qaoa
    Exp <- TRUE
}
#if(problem == 'qft')
#    state <- "sin"


source("experiment.R")
print(paste('Problem:',problem,'Graph:',graph,'Error:',Error,'error:',error,'e_mode:',e_mode,'tidx:',tidx))
print(paste('inputstate:',state,'Exp:',Exp,'depth:',depth,'RC:',RC,'extraCycles:',extraCycles,'samples:',samples))
state <- experiment(problem,graph,Error,error,e_mode,tidx,inputstate=state,Exp=Exp,
                    depth=depth,RC=RC,extraCycles=extraCycles,samples=samples)
