library(shiny)
library(ggplot2)
library(FITSio)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    # turn off annoying warnings
    options(warn=-1)
    
    observeEvent(c(input$thresh,input$strength,input$scale,input$filepath),{
        
        # check file exists
        filepath = input$filepath
        if(!file.exists(filepath)) filepath = paste0("data/",filepath)
        
        # load data and perform initial FFT
        galaxy = readFrameFromFITS(filepath)
        names(galaxy) = tolower(names(galaxy))
        f   = fft(galaxy$flux)
        
        ## initialize input list for debug purposes
        # input = list(strength=.02,thresh=100,scale="log",filepath="cB58_Lyman_break.fit")
        
        # define attenuation function using double exponential
        att = Vectorize(function(x){1/(exp(exp(input$strength*(x-input$thresh))))})
        
        # calculate weights using attenuation function
        wts   = att(1:length(f))
        
        # rescale according to weights
        f.att = (f * wts)/sum(wts)
        
        # invert FFT, then take real part
        galaxy$f.flt = Re(fft(f.att,inverse=T))
        
        # roughly rescale coefficients using linear regression on flux quantiles
        rescale_coefs = coef(lm(
            quantile(galaxy$flux,probs=c(seq(.15,.85,.05))) ~
                quantile(galaxy$f.flt,probs=c(seq(.15,.85,.05)))
        ))
        galaxy$f.flt = galaxy$f.flt * rescale_coefs[2] + rescale_coefs[1]
        
        # plot periodogram and attenuation curve
        output$freq <- renderPlot({
            ggplot(data.frame(mod=Mod(f)),aes(x=1:length(mod),y=mod)) +
                geom_line(aes(color="a")) + scale_y_log10(labels = function(x){signif(x,1)}) + 
                geom_smooth(aes(color="b"),se=F,method="loess",formula='y~x',span=.1) +
                stat_function(aes(color="c"),fun=function(x){att(x)+.1},size=1) + 
                scale_colour_manual(name="Legend",
                                    values=c("a"="black",
                                             "b"="red",
                                             "c"="blue"),
                                    labels=c("Original spectrum",
                                             "Approximate trend",
                                             "Attenuation curve")) +
                theme(legend.justification = c(1, 1), legend.position = c(.98, .98)) +
                labs(x="Frequency index",y="Relative intensity",title="Fourier transformed spectrum + attenuation curve")
        })
        
        # plot original spectrum
        output$cb58 <- renderPlot({
            switch(input$scale,
                   linear={
                       ggplot(galaxy,aes(x=loglam,y=flux)) + geom_line() + 
                           scale_y_continuous() + 
                           labs(x="Log(Wavelength)",y="Flux",
                                title="Original spectrum")
                   },
                   log={
                       ggplot(galaxy,aes(x=loglam,y=flux)) + geom_line() + 
                           scale_y_log10(limits=c(exp(min(log(galaxy$f.flt),na.rm=T)),max(galaxy$f.flt))) + 
                           labs(x="Log(Wavelength)",y="Log(Flux)",
                                title="Original spectrum")
                   })
            
            
        })
        
        # plot filtered spectrum
        output$cb58flt <- renderPlot({
            switch(input$scale,
                   linear={
                       ggplot(galaxy,aes(x=loglam,y=f.flt)) + geom_line() + 
                           labs(x="Log(Wavelength)",y="Flux",
                                title="Smoothed spectrum")
                   },
                   log={
                       ggplot(galaxy,aes(x=loglam,y=f.flt)) + geom_line() + 
                           scale_y_log10() + 
                           labs(x="Log(Wavelength)",y="Log(Flux)",
                                title="Smoothed spectrum")
                   })
        })
        
        # plot noise spectrum
        output$noise <- renderPlot({
            ggplot(galaxy,aes(x=loglam,y=flux-f.flt)) + geom_line() + 
                labs(x="Log(Wavelength)",y="Flux",
                     title="Noise (original-smoothed)")
        })
    })
})
