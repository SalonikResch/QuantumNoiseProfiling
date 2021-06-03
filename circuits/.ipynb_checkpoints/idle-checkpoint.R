#Write schedule, rather than make circuit
idle_ckt <- function(n,depth,schedule=FALSE){
	#circuit <- list()
	gates <- list()
    
    for(d in 1:depth)
        for(q in 1:n)
            gates <- c(gates,list(list('I',q-1,"")))

    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}

