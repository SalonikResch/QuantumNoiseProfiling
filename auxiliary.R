
source("simulate.R")
source("schedule.R")

get_circuit <- function(problem,graph="",schedule=FALSE,depth=10,RC=FALSE,extraCycles=FALSE){
	if(problem == 'qaoa'){
		source("circuits/qaoa.R")
		return(setupQAOA(graph,schedule=schedule,reorder=TRUE))
	}
	#Clifford+T QAOA
	if(problem == 'qaoaDecomposed'){
		source("circuits/qaoa.R")
		ckt <- setupQAOA(graph,schedule=schedule,reorder=TRUE,decompose=TRUE)
		if(RC)
			ckt <- randomizeCompile(ckt)
		if(extraCycles)
			ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
	}
	#Rotation Gate QFT
	if(problem == 'qft'){
		source("circuits/qft.R")
		return(qft_ckt(n=graph,schedule=schedule))
	}	
	#Clifford+T QFT
	if(problem == 'qftDecomposed'){
		source("circuits/qft.R")
		ckt <- qft_gridsynth_ckt(n=graph,schedule=schedule)
		if(RC)
			ckt <- randomizeCompile(ckt)
		if(extraCycles)
			ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
	}	
    #Adder
	if(problem == 'adder'){
		source("circuits/adder.R")
		ckt <- adder_ckt(graph,schedule=schedule)
		if(RC)
			ckt <- randomizeCompile(ckt)
		if(extraCycles)
			ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
	}
    if(problem == 'idle'){
        source("circuits/idle.R")
        ckt <- idle_ckt(n=graph,depth=depth,schedule=schedule)
        if(RC)
            ckt <- randomizeCompile(ckt)
        if(extraCycles)
            ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
    }
    if(problem == 'test'){
        source("circuits/testCircuits.R")
		return(test_ckt(idx=graph,schedule=schedule))
    }
    if(problem == 'adder_n4'){
        source("circuits/QASMBench.R")
        ckt <- adder_n4(schedule=schedule)
        if(RC)
            ckt <- randomizeCompile(ckt)
        if(extraCycles)
            ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
    }
    if(problem == 'linqaoa'){
        source("circuits/linqaoa.R")
        return(linqaoa(schedule=schedule,c=graph))
    }
    if(problem == 'GHZ'){
        source("circuits/GHZ.R")
        return(GHZ(nQubits=as.integer(graph),schedule=schedule))
    }
    if(problem == 'bv_n'){
        source("circuits/QASMBench.R")
        return(bv_n(nQubits=as.integer(graph),schedule=schedule))
    }
    if(problem == 'bv_nMOD'){
        source("circuits/QASMBench.R")
        return(bv_nMOD(nQubits=as.integer(graph),schedule=schedule))
    }
    if(problem == 'random'){
        source("circuits/random.R")
		ckt <- random_ckt(n=graph,depth=depth,schedule=schedule)
        if(RC)
            ckt <- randomizeCompile(ckt)
        if(extraCycles)
            ckt <- randomizeCompile(ckt,RC=FALSE)
		return(ckt)
    }
}

proceed <- function(circuit,e_mode,tidx,nPoints=100){
	#number of cycles is length of circuit list
	nCycles <- length(circuit)
    print(paste("Circuit has",nCycles,"cycles"))
	#Number of qubits is vertical dimension of circuit
	nQubits <- log(dim(circuit[[1]])[1],base=2)
	#number of gates is cycles times vertical dimension of circuit (# of gates in each cycle)
	nGates <- nCycles * nQubits
	
    #R setup indexes from 1, simulate.R converts to from 0 for QuantumOps
	if( (e_mode == 1 & tidx > nGates) | (e_mode == 2 & tidx >= nCycles) | (e_mode == 3 & tidx > nQubits) | (e_mode == 4 & tidx > nPoints) ){
		print("Parameters out of bounds")
		return(FALSE)
	}else
		return(TRUE)
}
    
get_scores <- function(problem,graph=""){
	#if(problem == 'qaoa'){
		source("circuits/qaoa.R")
		return(returnQAOAscores(graph))
	#}
}    



#Matrix: rows = qubits , columns = cycles
#Matrix for noise type, matrix for noise value
noise_setup <- function(circuit,e_mode,e_type,e_vals,tidx){

	nQubits <- log(dim(circuit[[1]])[1],base=2)	#get number of qubits
	nCycles <- length(circuit)

	e_typeMatrix <- matrix(0,nrow=nQubits,ncol=nCycles)
	e_valsMatrix <- matrix(0,nrow=nQubits,ncol=nCycles)

	#Convert names of noises to integer indices
	e_type[ which(e_type == 'AD') ] <- 1		#1 = Amplitude Damping
	e_type[ which(e_type == 'coherent') ] <- 2	#2 = coherent
	e_type[ which(e_type == 'pauli') ] <- 3		#3 = pauli
	e_type[ which(e_type == 'PD') ] <- 4		#4 = Phase Damping
    e_type[ which(e_type == 'coherentx') ] <- 5	#5 = coherentx
    e_type[ which(e_type == 'bothx') ] <- 6	#5 = coherentx
	e_type <- as.integer(e_type)
	

	#Create matrix of indices (organied by method of choice)
	if(e_mode == 0)
		idxMatrix <- matrix(0,nrow=nQubits,ncol=nCycles)
	#By gate
	if(e_mode == 1)
		idxMatrix <- matrix(1:(nQubits*nCycles),nrow=nQubits)
	#By cycle
	if(e_mode == 2)
		idxMatrix <- matrix(rep(1:nCycles,nQubits),nrow=nQubits,byrow=TRUE)
	#By qubit
	if(e_mode == 3)
		idxMatrix <- matrix(rep(1:nQubits,nCycles),nrow=nQubits)
	#By error level
	if(e_mode == 4)
		idxMatrix <- matrix(1,nrow=nQubits,ncol=nCycles)

	#tidx is target index for gates/cycles/qubits - but if changeError, then it was used for error level setting, now should be 1 (all)
	if(e_mode == 4)
		tidx <- 1
		
	#"Expand" if necessary
	if(length(e_mode == 1))
		e_mode <- rep(e_mode,length(tidx))
	if(length(e_type == 1))
		e_type <- rep(e_type,length(tidx))

	for(j in 1:length(tidx)){
		e_typeMatrix[ which(idxMatrix == tidx[j]) ] <- as.integer(e_type[j])
		e_valsMatrix[ which(idxMatrix == tidx[j]) ] <- as.numeric(e_vals[j])
	}

	list( e_typeMatrix , e_valsMatrix )
}



 
CSUR_write_data <- function(fidelity,problem,Error,error,graph,e_mode,tidx,tag="",Exp=FALSE,ExpVal=""){
	experiment_mode <- c('eachGate','eachCycle','eachQubit','changeError')[e_mode]

    #Fidelity
	if(problem == 'idle' || problem == 'random'){
		dir_fidelity <- paste("data/",problem,tag,"/fidelity/",experiment_mode,"/",graph,"/",Error,"/",error,sep="")
    		f <- paste(dir_fidelity,"/",paste(tidx,collapse=""),sep="")
	}else{
		dir_fidelity <- paste("data/",problem,tag,"/fidelity/",experiment_mode,"/",graph,"/",Error,sep="")
    		f <- paste(dir_fidelity,"/",paste(error,collapse=""),sep="")
	}
	dir.create(dir_fidelity,recursive=TRUE,showWarnings=FALSE)
    
	write(x=fidelity,file=f,ncolumns=1,append=TRUE)
	

	if(Exp){
        print(Exp)
		d <- paste("data/",problem,tag,"/ExpectationValue/",experiment_mode,"/",graph,"/",Error,sep="")
	        dir.create(d,recursive=TRUE,showWarnings=FALSE)
		f <- paste(d,"/",paste(error,collapse=""),sep="")
		
       		write(x=ExpVal,file=f,ncolumns=1,append=TRUE)
		
	}
}        
        
        
