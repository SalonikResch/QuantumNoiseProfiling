random_ckt <- function(n,depth,schedule=FALSE){
	#circuit <- list()
	gates <- list()
    
    easyGates <- c('I','X','Y','Z')
    easyIdx <- sample(1:4,size=ceiling(depth/2),replace=TRUE)
    
    for(d in 1:floor(depth/2)){
        #Easy Gates
        for(q in 1:n)
            gates <- c(gates,list(list(easyGates[easyIdx[d]],q-1,"")))
        #Hard gates (either one CX or many H)
        if(sample(1:2,size=1) == 1){ #CX
            gates <- c(gates,list(list('CX',sample(0:(n-1),size=2),"")))
        }else{ #H
            for(q in 1:n)
                gates <- c(gates,list(list('H',q-1,"")))
        }
    }
    #If odd # of cycles, add 1 extra easy cycle
    if(depth %% 2 == 1)
        for(q in 1:n)
            gates <- c(gates,list(list(easyGates[easyIdx[ceiling(depth/2)]],q-1,"")))       

    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=n,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=n,schedule(nQubits=n,gates=gates)))
}

