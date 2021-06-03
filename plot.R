
     
                    
#bottom <- 4.0; left <- 4.5; top <- 5.5; right <- 1.0;
#par(mar=c(bottom,left,top,right))                    
#legend("top",legend=legError,fill=colors,inset=c(0,-0.25),xpd=NA,cex=1.5,ncol=2,bty="n")
#legend("top",legend=c("Original","Randomly Compiled"),lty=c(1,3),lwd=3,inset=c(0.0,-0.12),xpd=NA,cex=1.3,ncol=2,bty="n")                    
                                  
CSUR_idleRandom <- function(){

    colors <- c('Blue','Green','Orange','Red')
    bottom <- 4.0; left <- 4.5; top <- 5.6; right <- 1.0;
    par(mar=c(bottom,left,top,right))
    es <- c(0,0.01,0.02,0.03)
    cs <- paste( round(seq(0,0.1,length.out=4),3),'*pi',sep='')
    ltext <- list(es,es,paste(es,cs,sep='/'),cs)             


    nQubits <- 4
    depth <- 2:35
    error <- c(0,0.01,0.02,0.03)
    Error <- c('AD','pauli','bothx','coherent')
    Error <- c('bothx')
    idleRandom <- c('idle','random')
    
    for(E in Error){ 
    for(IR in idleRandom){
    if(E == 'coherent')
        error <- c(0,0.103672557568463,0.210486707790516,0.314159265358979)
    folder <- paste("dataGraphs/",IR,sep='')
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    jpeg(paste(folder,'/',E,".jpg",sep=''))
    par(mar=c(bottom,left,top,right))
    plot(1:2,1:2,type='n',ylim=c(0,1),xlim=c(1,max(2*depth)),xlab='Cycles',ylab='Process Fidelity',cex.lab=1.25)
    pf <- rep(-1,length(depth))
    RCpf <- rep(-1,length(depth))
    pf95 <- rep(-1,length(depth))
    RCpf95 <- rep(-1,length(depth))
    for(e in error){
        for(d in depth){
            print(paste(E,IR,'depth',d,'error',e))
            idx <- which(depth == d)
            tpf <- rep(1,5)
            #100 indicates it was the "max error", the full value of error
            tpf <- as.numeric( readLines(paste('data/',IR,'/random/depth',d,'/extraCycles',
                                                         '/fidelity/changeError/',nQubits,'/',E,'/',e,'/100',sep='')) )
            
            trc <- as.numeric( readLines(paste('data/',IR,'/random/depth',d,'/RC',
                                                         '/fidelity/changeError/',nQubits,'/',E,'/',e,'/100',sep='')) )
            
            print(paste(' Normal:',length(tpf),'samples  RC:',length(trc),'samples'))
            pf[idx] <- mean(tpf)
            pf95[idx] <- 1.96*sd(tpf)/sqrt(length(tpf)) /2
            RCpf[idx]  <- mean(trc)
            RCpf95[idx] <- 1.96*sd(trc)/sqrt(length(trc)) /2
        }
        color <- colors[which(error == e)]
        lines(2*depth,  pf,lwd=4,lty=1,col=color)
        arrows(2*depth,pf+pf95,2*depth,pf-pf95,length=0.05,angle=90,code=3,col=color)
        lines(2*depth,RCpf,lwd=4,lty=3,col=color)
        arrows(2*depth,RCpf+RCpf95,2*depth,RCpf-RCpf95,length=0.05,angle=90,code=3,col=color)
    }
	print(ltext[[which(Error == E)]])
	legend("top",legend=ltext[[which(Error == E)]],fill=colors,inset=c(0,-0.22),xpd=NA,cex=1.5,ncol=2,bty="n")
	legend("top",legend=c("Original","Randomly Compiled"),lty=c(1,3),lwd=3,inset=c(0.0,-0.1),xpd=NA,cex=1.3,ncol=2,bty="n")
        print(pf)
    dev.off()
    }
    }
}
                    
                    
CSUR_changeError <- function(){
    EDIV <- 1 #Plot less than full error
    colors <- c('Black','Red','Green','Blue')
    bottom <- 4.0; left <- 4.5; top <- 5.6; right <- 1.0;
    
    Error <- c('AD','pauli','bothx','coherent')
    problem <- c('qft','qaoa')
    
    xlab <- list('Probability of Amplitude Damping','Pauli Error Rate','Pauli X Error Rate and X-Rotation Angle','Z-Rotation Angle')
    ylab <- list('Process Fidelity','Expectation Value')
    
    for(E in Error){ 
    for(P in problem){

    if(P == 'qft' || P == 'qftDecomposed'){
       metric <- 'fidelity'
       graph <- 4
       state <- 'random'
       ymax <- 1
    }
    if(P == 'qaoa' || P == 'qaoaDecomposed'){
       metric <- 'ExpectationValue'
       graph <- '3reg'
       state <- 'zero'
       ymax <- 12
    }

    #Setup output file
    folder <- paste("dataGraphs/",P,sep='')
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    jpeg(paste(folder,'/',E,".jpg",sep=''))
    par(mar=c(bottom,left,top,right))

    #Find all files (error levels) in directory and sort in increasing order
    folder <- paste('data/',P,'/',state,'/',metric,'/changeError/',graph,'/',E,sep='')
    error <- list.files(path=folder,full.names=FALSE,recursive=FALSE)
    error <- sort(as.numeric(error),decreasing=FALSE)

    print(paste("Folder ",folder,"has",length(error),"files"))
    #Read them all
    pf <- rep(-1,length(error))
    pf95 <- rep(-1,length(error))
    for(e in error){
            idx <- which(error == e)
            tpf <- rep(1,5)
            #100 indicates it was the "max error", the full value of error
            tpf <- as.numeric( readLines(paste(folder,e,sep='/')) )
            
            
            print(paste(P,E,'error',e,' Normal:',length(tpf),'samples '))
            pf[idx] <- mean(tpf)
            pf95[idx] <- 1.96*sd(tpf)/sqrt(length(tpf)) /2
    }
	print(pf)
    color <- colors[which(Error == E)]
        print(problem)
        print(P)
        print(which(problem == P))
        print(metric)
    plot(1:2,1:2,type='n',ylim=c(0,ymax),xlim=c(0,max(error)/EDIV),xaxt='n',xlab=xlab[[which(Error == E)]],ylab=ylab[[which(problem == P)]],cex.lab=1.5)
    lines(error,  pf,lwd=4,lty=1,col=color)
    #arrows(error,pf+pf95,error,pf-pf95,length=0.05,angle=90,code=3,col=color)
        
    ##Different axes for different plots
    serror <- round(error[floor(seq(from=1,to=length(error),length.out=7))]/EDIV,3)
    if(E == 'coherent')
        axis(1,at=serror,labels=paste(round(serror*3.33,3),'*pi',sep=''),cex=1.5,line=0,srt=45)
    else if(E == 'bothx')
        axis(1,at=serror,labels=paste(serror,paste(round(serror*3.33,3),'*pi',sep=''),sep='/'),cex=1.5,line=0)
    else
        axis(1,at=serror,labels=serror,cex=1.5,line=0)
    dev.off()
    }
    }
}
                    
CSUR_changeErrorRCandEC <- function(){
    EDIV <- 3 #Plot less than full error
    colors <- c('Black','Red','Green','Blue')
    bottom <- 4.0; left <- 4.5; top <- 5.6; right <- 1.0;

    Error <- c('AD','pauli','bothx','coherent')
    problem <- c('adder','qftDecomposed','qaoaDecomposed')
    
    xlab <- list('Probability of Amplitude Damping','Pauli Error Rate','Pauli X Error Rate and X-Rotation Angle','Z-Rotation Angle')
    ylab <- list('Process Fidelity','Process Fidelity','Expectation Value')
    
    for(E in Error){ 
    for(P in problem){


    if(P == 'qft' || P == 'qftDecomposed'){
       metric <- 'fidelity'
       graph <- 4
       state <- 'random'
       ymax <- 1
    }
    if(P == 'qaoa' || P == 'qaoaDecomposed'){
       metric <- 'ExpectationValue'
       graph <- '3reg'
       state <- 'zero'
       ymax <- 12
    }
    if( P == 'adder'){
        metric <- 'fidelity'
        graph <- 7
        state <- 'random'
        ymax <- 1
    }

    #Setup output file
    folder <- paste("dataGraphs/",P,sep='')
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    jpeg(paste(folder,'/',E,".jpg",sep=''))
    par(mar=c(bottom,left,top,right))


    #Find all files (error levels) in directory and sort in increasing order
    folder <- paste('data/',P,'/',state,'/extraCycles/',metric,'/changeError/',graph,'/',E,sep='')
    error <- list.files(path=folder,full.names=FALSE,recursive=FALSE)
    error <- sort(as.numeric(error),decreasing=FALSE)
    print(paste('There are',length(error),'files in folder',folder))

    #Read them all
    pf <- rep(-1,length(error))
    pf95 <- rep(-1,length(error))
    for(e in error){
            idx <- which(error == e)
            tpf <- as.numeric( readLines(paste(folder,e,sep='/')) )
            
            print(paste(P,E,'error',e,' Normal:',length(tpf),'samples '))
            pf[idx] <- mean(tpf)
            pf95[idx] <- 1.96*sd(tpf)/sqrt(length(tpf)) /2
    }

    #Again with RC
    folder <- paste('data/',P,'/',state,'/RC/',metric,'/changeError/',graph,'/',E,sep='')
    RCerror <- list.files(path=folder,full.names=FALSE,recursive=FALSE)
    RCerror <- sort(as.numeric(RCerror),decreasing=FALSE)
    print(paste('There are',length(RCerror),'files in folder',folder))

    #Read them all
    RCpf <- rep(-1,length(RCerror))
    RCpf95 <- rep(-1,length(RCerror))
    print(paste("Length of RCerror",length(RCerror),"length of RCpf",length(RCpf)))
        print(RCerror)
    for(e in RCerror){
            idx <- which(RCerror == e)
            #print(paste(e,'is the',idx,'element'))
            tpf <- as.numeric( readLines(paste(folder,e,sep='/')) )
            
            print(paste(P,E,'error',e,' RC:',length(tpf),'samples '))
            RCpf[idx] <- mean(tpf)
            RCpf95[idx] <- 1.96*sd(tpf)/sqrt(length(tpf)) /2
            print(paste("Min:",min(tpf),"Mean:",mean(tpf),"Max:",max(tpf),"sd:",sd(tpf)))
    }
        
    #Plot
    color <- colors[which(Error == E)]
    plot(1:2,1:2,type='n',ylim=c(0,ymax),xlim=c(0,max(error)/EDIV),xaxt='n',xlab=xlab[[which(Error == E)]],ylab=ylab[[which(problem == P)]],cex.lab=1.5)
    lines(error,  pf,lwd=4,lty=1,col=color)
    #arrows(error,pf+pf95,error,pf-pf95,length=0.05,angle=90,code=3,col=color)
        print(paste(length(RCerror),length(RCpf)))
    lines(RCerror,  RCpf,lwd=4,lty=3,col=color)
        print("done")
    #arrows(RCerror,RCpf+RCpf95,error,RCpf-RCpf95,length=0.05,angle=90,code=3,col=color)
    legend("top",legend=c("Original","Randomly Compiled"),lty=c(1,3),lwd=3,inset=c(0.0,-0.1),xpd=NA,cex=1.3,ncol=2,bty="n")
        
    ##Different axes for different plots
    serror <- round(error[floor(seq(from=1,to=length(error),length.out=7))]/EDIV,3)
    if(E == 'coherent')
        axis(1,at=serror,labels=paste(round(serror*3.33,3),'*pi',sep=''),cex=1.5,line=0)
    else if(E == 'bothx')
        axis(1,at=serror,labels=paste(serror,paste(round(serror*3.33,3),'*pi',sep=''),sep='/'),cex=1.5,line=0)
    else
        axis(1,at=serror,labels=serror,cex=1.5,line=0)
    dev.off()
    }
    }
}
