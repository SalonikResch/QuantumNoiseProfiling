addGate <- function(gates,gate,qubits,arguments=""){
    return( c(gates, list( list(gate,qubits,arguments) )) )
}


#Each gate has a list associated with it
#[[ Name , [qubits] , parameter ]]

#input is a list of these lists
#output is same, but organized by cycle

schedule <- function(nQubits,gates){

	#vector of when qubits are available next
	available <- rep(1,nQubits)

	#schedule is a list, with a list for each cycle
	schedule <- list()

	for(g in gates){
		#soonest that all qubits are available
		soonest <- max(available[g[[2]]+1])			

		#schedule it
		if(length(schedule) >= soonest && length(schedule[[soonest]]) > 0)
			schedule[[soonest]] <- c(schedule[[soonest]],list(g))		
		else
			schedule[[soonest]] <- list(g)

		#then indicate qubits are busy until after
		available[g[[2]]+1] <- soonest+1						
		#print(paste(g[[1]],"on cycle",soonest))
	}
	#printSchedule(nQubits,schedule,plot=FALSE)
	return(schedule)	
}


schedule2circuit <- function(nQubits,schedule){
	circuit <- list()
    print("nQubits")
    print(nQubits)

	#Iterate over cycles
	for(c in 1:length(schedule)){
		m <- many(gate=I(),n=nQubits)	#Start with Identity for each cycle

		#Build matrix for cycle by multiply by matrix for each gate
		for(g in schedule[[c]]){
			#Finite number of gates, so just check for each
			if(g[[1]] == 'Rz')
				m <- single(gate=Rz(theta=g[[3]]),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'Rx')
				m <- single(gate=Rx(theta=g[[3]]),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'CX')
				m <- controlled(gate=X(),n=nQubits,cQubits=g[[2]][1],tQubit=g[[2]][2]) %*% m
			else if(g[[1]] == 'CRz')
				m <- controlled(gate=Rz(theta=g[[3]]),n=nQubits,cQubits=g[[2]][1],tQubit=g[[2]][2]) %*% m
			else if(g[[1]] == 'u1')
				m <- single(gate=u1(theta=g[[3]]),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'u2')
				m <- single(gate=u2(phi=g[[3]][1],lambda=g[[3]][2]),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'u3')
				m <- single(gate=u3(theta=g[[3]][1],phi=g[[3]][2],lambda=g[[3]][3]),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'Cu1')
				m <- controlled(gate=u1(theta=g[[3]]),n=nQubits,cQubits=g[[2]][1],tQubit=g[[2]][2]) %*% m

			else if(g[[1]] == 'H')
				m <- single(gate=H(),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'X')
				m <- single(gate=X(),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'Y')
				m <- single(gate=Y(),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'Z')
				m <- single(gate=Z(),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'S')
				m <- single(gate=S(),n=nQubits,t=g[[2]]) %*% m
			else if(g[[1]] == 'T')
				m <- single(gate=T(),n=nQubits,t=g[[2]]) %*% m
            else if(g[[1]] == "T'")
				m <- single(gate=adjoint(T()),n=nQubits,t=g[[2]]) %*% m
            else if(g[[1]] == 'I')
				m <- single(gate=I(),n=nQubits,t=g[[2]]) %*% m  #Does nothing, just for consistency
			else
				print(paste("Did not recognize gate:",g[[1]]))
		}
		circuit <- c(circuit, list(m) )
	}
	return(circuit)
}


printSchedule <- function(nQubits,schedule,plot=FALSE){
	colors <- c("Black","Red")
    if(plot)
        plot(1:length(schedule),1:length(schedule),type="n",ylim=c(0,nQubits-1))
	for(c in 1:length(schedule)){
		print(paste("Cycle ",c,":",sep=""))
		for(g in schedule[[c]]){
			print(paste( "Qubit(s)", paste(g[[2]],collapse=' '),g[[1]],g[[3]]))
            if(plot){
                for(j in 1:length(g[[2]]))
                    text(c,g[[2]][j],g[[1]],col=colors[j])
                if(length(g[[2]]) > 1)
                    lines(c(c,c),c(g[[2]][1],g[[2]][2]))
            }
		}
		print("")
	}
}


        
##Randomized Compiling (Matrix version, cannot be used with schedule as explicit gates)      
#This assumes that all gates in circuit are compatible with RC (Clifford+T), does not check
        
#Interleave circuit with RC cycles (or idle gates)
randomizeCompile <- function(circuit,RC=TRUE){
    #Get # of qubits
    nQubits <- log(dim(circuit[[1]])[1],base=2)
    #Get # of cycles
    nCycles <- length(circuit)
    
    if(!RC){
        m <- I()
        for(k in 2:nQubits)
            m <- tensor(m,I())
        newCircuit <- list()
        for(j in 1:length(circuit)){
            newCircuit <- c(newCircuit,list(m))
            newCircuit <- c(newCircuit,list(circuit[[j]]))
        }
        newCircuit <- c(newCircuit,list(m))
        return(newCircuit)
    }
    
    #Generate the randomize cycles
    PauliGates <- list(I(),X(),Y(),Z())  #Pauli Gates to choose from
    Twirling <- list()                   #Twirling gates (cycles)
    for(j in 1:nCycles){
        p <- sample(1:4,size=nQubits,replace=TRUE) #pick twirl gates
        m <- PauliGates[[p[1]]]   #Build Twirling matrix
        for(k in 2:nQubits)
            m <- tensor(m,PauliGates[[p[k]]])
        #Add it to list
        Twirling <- c(Twirling, list(m) )
    }
    #Find the corrective cycles
    Tc <- list() #Corrective gates
    for(j in 1:nCycles)
        Tc <- c(Tc, list( circuit[[j]] %*% adjoint(Twirling[[j]]) %*% adjoint(circuit[[j]]) ))

    #Create new circuit
    #Start
    RCcircuit <- list(Twirling[[1]])                 #Start with Twirling
    RCcircuit <- c(RCcircuit, list(circuit[[1]]) )   #Then first original cycle ("Hard")
    #Middle    
    for(j in 2:nCycles){
        #Undo previous Twirling, add new twirling
        m <- Twirling[[j]] %*% Tc[[j-1]]
        RCcircuit <- c(RCcircuit,list(m))
        #Do next "Hard" cycles
        RCcircuit <- c(RCcircuit,list(circuit[[j]]))
    }
    #End
    RCcircuit <- c(RCcircuit,list(Tc[[nCycles]]))
    return(RCcircuit)
}

        
circuit2matrix <- function(circuit){
    m <- circuit[[1]]
    for(j in 2:length(circuit))
        m <- circuit[[j]] %*% m
    return(m)
}
        
        
        
        
        
        
        
        
        
        
