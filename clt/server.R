library(shiny)
library(sn)

# Define server logic required to draw a histogram
shinyServer(function(input,output) {
	
	observeEvent(input$dimension,{
		
		output$distPlot <- renderPlot(height=input$dimension[1],width=input$dimension[2],{
			
			size = c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,90,100,150,200,250,300,350,400,600,800,1000,2000,5000,10000)[input$size+1]
			type = input$dist
			S    = 10000
			means = rep(NA,S)
			
			switch(type,
						 
						 
						 normal={
						 	
						 	for(i in 1:S){
						 		samp = rnorm(size,0,1)
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(-3,3,.01)
						 	y = dnorm(x,0,1)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Normal distribution")*";   "*mu*" = 0 , "*sigma*" = 1"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 bimodal={
						 	
						 	for(i in 1:S){
						 		n1 = rbinom(1,size,.5) ; n2 = size-n1
						 		samp = c(rnorm(n1,-1,.4),rnorm(n2,1,.4))
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(-3,3,.01)
						 	y = (dnorm(x,-1,.4)+dnorm(x,1,.4))/2
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Bimodal")*";   N(-1,0.4) + N(1,0.4)"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 skewed={
						 	
						 	a = 10
						 	w = sqrt(pi/2*(1+1/a^2))
						 	
						 	for(i in 1:S){
						 		samp = rsn(size,-1,w,a)
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(-3,3,.01)
						 	y = dsn(x,-1,w,a)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Skewed normal")*";   "*xi*" = -1 , "*alpha*" = 10 , "*omega*" = "*sqrt(pi*(1+alpha^-2)/2)),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 uniform={
						 	
						 	for(i in 1:S){
						 		samp = runif(size,0,1)
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(-.5,1.5,.002)
						 	y = dunif(x,0,1)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Uniform")*";   a = 0 , b = 1"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 weibull={
						 	
						 	for(i in 1:S){
						 		samp = rweibull(size,0.5)
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(0,4,.002)
						 	y = dweibull(x,0.5)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Weibull")*";   "*lambda*" = 1 , k = 0.5"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 beta={
						 	
						 	for(i in 1:S){
						 		samp = rbeta(size,0.5,0.5)
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(0,1,.002)
						 	y = dbeta(x,0.5,0.5)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Beta")*";   "*alpha*" = "*beta*" = 0.5"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 },
						 
						 
						 
						 sawtooth={
						 	
						 	for(i in 1:S){
						 		samp = runif(size,0,0.5)+runif(size,0,0.5)
						 		samp = c(samp[samp<0.5],1.5-samp[samp>0.5])
						 		means[i] = mean(samp)
						 	}
						 	
						 	x = seq(0,1.002,.002)
						 	vy = Vectorize({function(x) if(x>0 && x<=0.5){4*x} else if(x>0.5 && x<=1){4*x-2} else 0})
						 	y = vy(x)
						 	ymax = c(y,hist(samp,plot=FALSE)$density)
						 	ymax = max(ymax[is.finite(ymax)])
						 	
						 	par(mfrow=c(3,1),cex=1,bty="l",oma=c(2,2,-.1,2)+0.1,mar=c(2.5,4,2.5,0)+0.1)
						 	plot(x,y,main=expression(bold("Sawtooth")*";   "*mu*" = 7/12 "%~~%"0.58"),ylab="Density",xlab="",type='l',lwd=2)
						 	hist(samp,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample (",size," observations)",sep=""),ylim=c(0,ymax),xlab="")
						 	lines(x,y,lty=2,lw=1.5,col="red")
						 	hist(means,freq=FALSE,xlim=c(min(x),max(x)),main=paste("Distribution of sample mean (",S," samples)",sep=""),xlab="")
						 }
			)
		})			 
	})
})