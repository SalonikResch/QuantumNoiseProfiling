

simulate <- function(circuit,state="",e_mode=0,e_type="",e_vals=""){
	
	for(c in 1:length(circuit)){
		#Perform gate operation
		state <- circuit[[c]] %*% state %*% adjoint(circuit[[c]])

		#Apply the noise  
		if(e_mode != 0){
			ADtargets <- which(e_type[,c] == 1) #Amplitude Damping
			CHtargets <- which(e_type[,c] == 2) #Coherent (Z)
			PNtargets <- which(e_type[,c] == 3) #Pauli
			PDtargets <- which(e_type[,c] == 4) #Phase Damping
            XCtargets <- which(e_type[,c] == 5) #Coherent (X)
            BXtargets <- which(e_type[,c] == 6) #Combine multiple for BothX
			if(length(ADtargets) > 0)
				state <- AmplitudeDamping(	p=state,	Pad=e_vals[ADtargets,c],	t=ADtargets-1)
			if(length(CHtargets) > 0)
				state <- CoherentNoise(		p=state,	theta=e_vals[CHtargets,c],	t=CHtargets-1)
            #Pauli noise currently not as general as others, same magnitude must be applied to each of the targets
			if(length(PNtargets) > 0)
			#	state <- PauliNoise(		p=state,	e=e_vals[PNtargets,c],		t=PNtargets)		PauliNoise needs some work
                state <- PauliNoise(		p=state,	e=e_vals[PNtargets[1],c],	t=PNtargets)
			if(length(PDtargets) > 0)
				state <- PhaseDamping(		p=state,	Ppd=e_vals[PDtargets,c],	t=PDtargets-1)
                
            if(length(XCtargets) > 0)
                state <- CoherentNoise(	p=	state,		theta=e_vals[XCtargets,c],	t=XCtargets-1,type='X') 
                
            # - changing error rate for coherent, to match with PauliX
            if(length(BXtargets) > 0){
                state <- CoherentNoise(		p=state,	theta=e_vals[BXtargets,c]*3.33*pi,	t=BXtargets-1)
                state <- PauliNoise(		p=state,	ex=e_vals[BXtargets[1],c],      t=BXtargets, ey=0,ez=0)
            }
		}

	}
	return(state)
}
















