reorder_edges <- function(edges){
    newEdges <- matrix(-1,nrow=dim(edges)[1],ncol=2)
    newEdges[1,] <- edges[1,]
    dirty <- edges[1,]
    completed <- 1
    #Iterate until all completed
    while(length(completed) < dim(edges)[1]){
        repeat{
            #Find rows (edges) which do not contain dirty qubits
            avail <- which( apply(X=edges,FUN=function(x)!any(x%in%dirty),MARGIN=1) )
            #Remove those which are already completed
            avail <- avail[ ! avail %in% completed ]
            #Check if none available
            if(length(avail) == 0)
                break
            #Add edges to new edges
            newEdges[length(completed)+1,] <- edges[avail[1],]
            #Add row index to completed
            completed <- c(completed,avail[1])
            #Add both qubits to dirty list
            dirty <- c(dirty,edges[avail[1],])
        }
        dirty <- vector()
    }
    return(newEdges)
}



#Write schedule, rather than make circuit
qaoa_ckt <- function(n,p,beta,gamma,edges,schedule=FALSE,reorder=FALSE){
    if(reorder)
        edges <- reorder_edges(edges)
	#circuit <- list()
	gates <- list()

	#Initial Hadamards
	#circuit <- c(circuit,list( many(gate=H(),n=n)))
	for(j in 1:n)
		gates <- c(gates,list(list('H',j-1,"")))

	#For the number of specified iterations
	for( P in 1:p){
		#Phase 1
		#For each edge
		for(j in 1:dim(edges)[1]){
			#First CNOT
			#circuit <- c(circuit,list( controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) ))
			gates <- c(gates,list(list('CX',edges[j,],"")))

			#Rz on target qubit
			#circuit <- c(circuit,list( single(gate=Rz(theta=gamma[P]),n=n,t=edges[j,2]) ))
			gates <- c(gates,list(list('Rz',edges[j,2],gamma[P])))

			#Second CNOT
			#circuit <- c(circuit,list( controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) ))
			gates <- c(gates,list(list('CX',edges[j,],"")))
		}

		#Phase 2
		#Fine x-rotations ARE possible
		#circuit <- c(circuit,list( many(gate=Rx(theta=2*beta[P]),n=n) ))
		for(j in 1:n)
			gates <- c(gates,list(list('Rx',j-1,2*beta[P])))
	}
    
    for(i in 1:length(gates))
        print(paste(gates[[i]]))
    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}

#condensed version of circuit (if point is to evaluate output, not simulate)
qaoa_matrix <- function(n,p,beta,gamma,edges,reorder=FALSE){
    if(reorder)
        edges <- reorder_edges(edges)
	m <- diag(2^n)
	m <- many(gate=H(),n=n) %*% m 	#Initial Hadamards
	for( P in 1:p){			#For the number of specified iterations
		for(j in 1:dim(edges)[1]){
			m <- controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) %*% m 	#First CNOT
			m <- single(gate=Rz(theta=gamma[P]),n=n,t=edges[j,2]) %*% m					#Rz on target qubit
			m <- controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) %*% m	#Second CNOT
		}
		m <- many(gate=Rx(theta=2*beta[P]),n=n) %*% m									#Fine x-rotations ARE possible
	}
	return(m)
}

#IBM qaoa	
qaoa_IBM_ckt <- function(n,beta,gamma,edges){
    edges <- n - edges - 1
    gates <- list()
	#Initial Hadamards
	for(j in 1:n)
		gates <- c(gates,list(list('H',j-1,"")))
	for(j in 1:dim(edges)[1]){
		#Controlled u1 (Rz?)
		gates <- c(gates,list(list('Cu1',edges[j,],-2*gamma)))
		#u1 (Rz?) on both qubits
		gates <- c(gates,list(list('u1',edges[j,1],gamma)))
		gates <- c(gates,list(list('u1',edges[j,2],gamma)))
	}
	#Rxs
	for(j in 1:n)
		gates <- c(gates,list(list('Rx',j-1,2*beta)))
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}

#qaoa_IBM_matrix <- function(n,beta,gamma,edges){
#	m <- diag(2^n)
#	m <- many(gate=H(),n=n) %*% m 	#Initial Hadamards
#	for(j in 1:dim(edges)[1]){
#		m <- controlled(gate=u1(theta=-2*gamma),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) %*% m 	#First CNOT#
#		m <- single(gate=u1(theta=gamma),n=n,t=edges[j,1]) %*% m					#Rz on control qubit
#		m <- single(gate=u1(theta=gamma),n=n,t=edges[j,2]) %*% m					#Rz on target qubit
#	}
#	m <- many(gate=Rx(theta=2*beta),n=n) %*% m									#Fine x-rotations ARE possible
#	return(m)
#}	

qaoa_IBM_matrix <- function(n,beta,gamma,edges){
    circuit <- qaoa_IBM_ckt(n,beta,gamma,edges)
    g <- circuit[[1]]
    for(j in 2:length(circuit))
        g <- circuit[[j]] %*% g
    return(g)
}

#Clifford+T QAOA
getDec <- function(angle){
    if(angle < 0)
        angle <- paste('\"(',angle,')\"',sep='')
    number <- sample(1:1000000000,size=1)
    sFile <- paste("sequence",number,sep="")
    cmd <- paste("./gridsynth ",angle," -b ",4," > ",sFile,sep="")	#store in sequence
	print(cmd)
    system(cmd)								#external call
    gates <- readLines(sFile)			#read sequence
    file.remove(sFile)
    gates <- rev(strsplit(gates,"")[[1]])
    gates <- gates[ gates %in% c('H','T','S','X') ]
    return(gates)
}
qaoa_ckt_decomposed <- function(n,p,beta,gamma,edges,schedule=FALSE,reorder=FALSE){
	#There is a file associated with each decomposed qaoa circuit
	cFolder <- "circuits/precompiled/qaoa"
	cFile <- paste(cFolder,paste(edges,collapse=''),sep='/')

	#if it exists, just read it in
	if(file.exists(cFile)){
		print("Reading circuit file")
		gates <- list()
		Lines <- readLines(cFile)
		for(L in Lines){
			L <- unlist(strsplit(L,split=' '))
			gate <- L[1]
			qubits <- as.integer(unlist(strsplit(L[2],split=',')))
			if(length(L) > 2)
				args <- as.numeric(unlist(strsplit(L[3],split=',')))
			else 
				args <- ''
			gates <- c(gates,list(list(gate,qubits,args)))
		}
	}else{
	#else, create it and write it
		print("Creating and writing circuit")	
	        dir.create(cFolder,recursive=TRUE,showWarnings=FALSE)
		if(reorder)
			edges <- reorder_edges(edges)
		#circuit <- list()
		gates <- list()

		#Initial Hadamards
		#circuit <- c(circuit,list( many(gate=H(),n=n)))
		for(j in 1:n)
		gates <- c(gates,list(list('H',j-1,"")))

		#For the number of specified iterations
		for( P in 1:p){
			#Phase 1
			#For each edge
			for(j in 1:dim(edges)[1]){
				#First CNOT
				#circuit <- c(circuit,list( controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) ))
				gates <- c(gates,list(list('CX',edges[j,],"")))
	
				#Rz on target qubit
				lg <- getDec(gamma[P]);   if(length(lg) < 2) print("Error: Not enough precision on gridsynth");
				for(i in 1:length(lg))
					gates <- c(gates,list(list(lg[i],edges[j,2],"")))		

				#Second CNOT
				#circuit <- c(circuit,list( controlled(gate=X(),n=n,cQubits=edges[j,1],tQubit=edges[j,2]) ))
				gates <- c(gates,list(list('CX',edges[j,],"")))
			}

			#Phase 2
			for(j in 1:n){
				gates <- c(gates,list(list('H',j-1,'')))
				lg <- getDec(2*beta[P]);   if(length(lg) < 2) print("Error: Not enough precision on gridsynth");
				for(i in 1:length(lg))
					gates <- c(gates,list(list(lg[i],j-1,"")))		
				gates <- c(gates,list(list('H',j-1,'')))
			}
		}
		#Write the gates
		for(g in gates){
			write(x=paste(g[[1]],paste(g[[2]],collapse=','),paste(g[[3]],collapse=','),sep=' '),file=cFile,append=TRUE)
		}
	}
    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}


evaluateAllScores <- function(n,edges){
	edges <- edges + 1
	scores <- rep(0,2^n)	#For all possible bit strings (assignments)
	for(j in 1:(2^n)){
		b <- convert_dec2bin(j-1,n)	#R from 1		#get binary strings
		for(k in 1:(dim(edges)[1])){				#See if each edge generates a point
			if( b[edges[k,1]] != b[edges[k,2]] )
				scores[j] <- scores[j] + 1
		}
	}
	scores
}



#Read input functions
readN <- function(name){
	f <- paste("circuits/graphs/",name,"/n",sep="")
	n <- as.integer(readLines(f))
	n
}
readP <- function(name){
	f <- paste("circuits/graphs/",name,"/p",sep="")
	p <- as.integer(readLines(f))
	p	
}
readBeta <- function(name){
	f <- paste("circuits/graphs/",name,"/beta",sep="")
	beta <- as.numeric(readLines(f))
	beta
}
readGamma <- function(name){
	f <- paste("circuits/graphs/",name,"/gamma",sep="")
	gamma <- as.numeric(readLines(f))
	gamma
}
readEdges <- function(name){
	f <- paste("circuits/graphs/",name,"/edges",sep="")
	d <- matrix(as.integer(readLines(f)),ncol=2)
	d
}

setupQAOA <- function(graph,schedule=FALSE,reorder=FALSE,decomposed=FALSE){
	n <- readN(graph)
	p <- readP(graph)
	beta <- readBeta(graph)
	gamma <- readGamma(graph)
	edges <- readEdges(graph)

	if(decomposed)
		return(qaoa_ckt_decomposed(n=n,p=p,beta=beta,gamma=gamma,edges=edges,schedule=schedule,reorder=reorder))
	else
		return(qaoa_ckt(n=n,p=p,beta=beta,gamma=gamma,edges=edges,schedule=schedule,reorder=reorder))
}
    
returnQAOAscores <- function(graph){
	n <- readN(graph)
	edges <- readEdges(graph)

	return(evaluateAllScores(n=n,edges=edges))
}
