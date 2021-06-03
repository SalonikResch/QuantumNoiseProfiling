#Write schedule, rather than make circuit
qft_ckt <- function(n,schedule=FALSE){
	#circuit <- list()
	gates <- list()

	for(j in 0:(n-1)){
		gates <- c(gates,list(list('H',j,"")))	#H on each qubit
		if(j < n-1)									#Only do if not last qubit
		for(k in (j+1):(n-1)){
			angle <- pi / 2^(k-j)
			gates <- c(gates,list(list('u1',j, angle/2)))	#Rz on target by  angle/2
			gates <- c(gates,list(list('CX',c(k,j),"")))	#CNOT from control to target
			gates <- c(gates,list(list('u1',j,-angle/2)))	#Rz on target by -angle/2
			gates <- c(gates,list(list('CX',c(k,j),"")))	#CNOT from control to target
			#gates <- c(gates,list(list('Cu1',c(k,j),angle)))
		}
	}
    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}

qft_matrix <- function(n){
	circuit <- qft_ckt(n)
	g <- circuit[[1]]
	for(j in 2:length(circuit))
		g <- circuit[[j]] %*% g
	return(g)
}

qft_m <- function(n){
	m <- diag(2^n)
	for(j in 0:(n-1)){
		m <- single(H(),n=n,t=j) %*% m
		if(j < n-1)									#Only do if not last qubit
		for(k in (j+1):(n-1)){
			angle <- pi / 2^(k-j)
			m <- controlled(R(angle),n=n,cQubits=k,tQubit=j) %*% m
		}
	}
	return(m)
}

 
getDec <- function(angle){
    if(angle < 0)
        angle <- paste('\"(',angle,')\"',sep='')
    number <- sample(1:1000000000,size=1)
    sFile <- paste("sequence",number,sep="")
    cmd <- paste("./gridsynth ",angle," -b ",8," > ",sFile,sep="")	#store in sequence
	print(cmd)
    system(cmd)								#external call
    gates <- readLines(sFile)			#read sequence
    file.remove(sFile)
    gates <- rev(strsplit(gates,"")[[1]])
    gates <- gates[ gates %in% c('H','T','S','X') ]
    return(gates)
}
    
qft_gridsynth_ckt <- function(n,schedule=FALSE){
	#circuit <- list()
	gates <- list()

	for(j in 0:(n-1)){
		gates <- c(gates,list(list('H',j,"")))	#H on each qubit
		if(j < n-1)									#Only do if not last qubit
		for(k in (j+1):(n-1)){
			angle <- pi / 2^(k-j)
			#Rz on target by  angle/2
    	        lg <- getDec(angle/2)
		if(length(lg) < 2)
			print("Error: Not enough precision on gridsynth")
        	for(i in 1:length(lg))
                	gates <- c(gates,list(list(lg[i],j,"")))
		gates <- c(gates,list(list('CX',c(k,j),"")))	#CNOT from control to target
		#Rz on target by -angle/2
            	lg <- getDec(-angle/2)
            	for(i in 1:length(lg))
                	gates <- c(gates,list(list(lg[i],j,"")))    
		gates <- c(gates,list(list('CX',c(k,j),"")))	#CNOT from control to target
		}
	}
    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}    

    
