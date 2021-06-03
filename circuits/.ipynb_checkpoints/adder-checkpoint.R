Toffoli <- function(a,b,c){
    gates <- list()
    gates <- c(gates,list(list('H',c,'')))
    gates <- c(gates,list(list('CX',c(b,c),'')))
    gates <- c(gates,list(list("T'",c,'')))
    gates <- c(gates,list(list('CX',c(a,c),'')))
    gates <- c(gates,list(list('T',c,'')))
    gates <- c(gates,list(list('CX',c(b,c),'')))
    gates <- c(gates,list(list("T'",c,'')))
    gates <- c(gates,list(list('CX',c(a,c),'')))
    gates <- c(gates,list(list('T',b,'')))
    gates <- c(gates,list(list('T',c,'')))
    gates <- c(gates,list(list('CX',c(a,b),'')))
    gates <- c(gates,list(list('H',c,'')))
    gates <- c(gates,list(list('T',a,'')))
    gates <- c(gates,list(list("T'",b,'')))
    gates <- c(gates,list(list('CX',c(a,b),'')))
    gates
}

FullAdd <- function(cin,a,b,cout){
    gates <- list()
    gates <- c(gates,Toffoli(a=a,b=b,c=cout))
    gates <- c(gates,list(list('CX',c(a,b),'')))
    gates <- c(gates,Toffoli(a=cin,b=b,c=cout))
    gates <- c(gates,list(list('CX',c(cin,b),'')))
    gates
}


adder_ckt <- function(n,schedule=FALSE){
    nQubits <- n   #
    n <- (n-1)/3
    
	#circuit <- list()
	gates <- list()

	for(j in 1:n){
        idx <- 3*(j-1)
		gates <- c(gates,FullAdd(cin=idx,a=idx+1,b=idx+2,cout=idx+3))
	}
    #If just want the schedule (for graphing purposes)
	if(schedule)
        return(schedule(nQubits=nQubits,gates=gates))
    #Normally, get schedule and make a circuit
	return(schedule2circuit(nQubits=nQubits,schedule(nQubits=nQubits,gates=gates)))
}
