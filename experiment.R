
#All required functions
source("noiseAnalysis/auxiliary.R")

experiment <- function(problem,graph,Error,error,e_mode,tidx,inputstate="",NoiseInjection="Default",Exp=FALSE,
                       depth=-1,RC=FALSE,extraCycles=FALSE,samples=1){
    
    
	#Get circuit
	circuit <- get_circuit(problem=problem,graph=graph,depth=depth,RC=RC,extraCycles=extraCycles)
    
	#check if simulation should proceed
	if(!proceed(circuit=circuit,e_mode=e_mode,tidx=tidx))
		return(0)
        
    #Create Input State
    nQubits = log(dim(circuit[[1]])[1],base=2)	#get number of qubits
 
##
#Do multiple samples here
F <- rep(NA,samples); HF <- rep(NA,samples); counts <- matrix(NA,nrow=samples,ncol=2^nQubits); 
PST <- rep(NA,samples); ExpVal <- rep(NA,samples);
for(s in 1:samples){	        
##        
    
	if(inputstate == "" | inputstate == "zero")
		state <- convert_ket2DM(intket(x=0,n=nQubits))
	else if(inputstate == 'h')
		state <- convert_ket2DM(do.call(ket,as.list(rep(1,2^nQubits))))
    else if(inputstate == 'sin')
        state <- convert_ket2DM(do.call(ket,as.list( sin((1:2^nQubits)/5) )))
    else if(inputstate == 'random')
        state <- convert_ket2DM(ranket(nQubits))

    #Get the DM of the ideal operation
    ideal_state <- simulate(circuit=circuit,e_mode=0,state=state)

    #If standard mode of noise injection (by gate, cycle, or qubit)
    if(NoiseInjection == "Default"){
        #Create Noise matrices
        L <- noise_setup(circuit=circuit,e_mode=e_mode,e_type=Error,e_vals=error,tidx=tidx)
        e_type <- L[[1]]
        e_vals <- L[[2]]

        #Run circuit with noise impacts
        noisy_state <- simulate(circuit=circuit,e_mode=e_mode,e_type=e_type,e_vals=e_vals,state=state)
    }else if(NoiseInjection == "Rotation"){
        #Noise is inserted by over/under rotations for each gate
        noisy_circuit <- get_circuit(problem=problem,graph=graph,tidx=tidx,error=error)
        noisy_state <- simulate(circuit=noisy_circuit,e_mode=0,state=state)
    }
	F[s] <- fidelity(ideal_state,noisy_state)    
	HF[s] <- hellinger_fidelity(ideal_state,noisy_state) 
    counts[s,] <- Re(diag(noisy_state))
    PST[s] <- counts[s,which(counts[s,] == max(counts[s,]))][1]
    ExpVal[s] <- 0
    if(Exp){
        ExpVal[s] <- sum( counts[s,] * get_scores(problem=problem,graph=graph) )
        print(paste("Exp:",ExpVal[s]))
    }  

    print(paste("fidelity:",F[s],"hellinger_fidelity:",HF[s]))
    print(paste("Finished sample",s,"of",samples))
}      
        
    #Create a tag to add depth,RC,extraCycles specifiers (if they apply)
    tag <- paste( c('','/')[as.integer(any(c(depth>0,RC,extraCycles,samples >=0 )))+1],
		  c(inputstate),
                  c('',paste("/depth",depth,sep=''))[as.integer(depth > 0)+1],
                  c('','/RC')[as.integer(RC)+1],
                  c('','/extraCycles')[as.integer(extraCycles)+1],
                  sep='')
   if(FALSE){
	write_data(counts=counts,fidelity=F,hellinger_fidelity=HF,PST=PST,ExpVal=ExpVal,Exp=Exp,
               fidelitySD=FSD,hellinger_fidelitySD=HFSD,PSTSD=PSTSD,ExpValSD=ExpValSD,
               problem=problem,Error=Error,error=error,graph=graph,e_mode=e_mode,tidx=tidx,tag=tag)
    }
    CSUR_write_data(fidelity=F,problem=problem,Error=Error,error=error,graph=graph,e_mode=e_mode,tidx=tidx,tag=tag,Exp=Exp,ExpVal=ExpVal)
}




