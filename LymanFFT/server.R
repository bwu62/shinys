library(shiny)
library(ggplot2)
library(FITSio)


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    # turn off annoying warnings
    options(warn=-1)
    
    # declare reactive values
    steps    = reactiveVal(NA)
    spectrum = reactiveVal(NA)
    unsaved  = reactiveVal(NA)
    
    # run analysis on file
    observeEvent(c(input$filepath,input$thresh,input$strength),{
        
        # find file
        filepath = input$filepath
        if(!file.exists(filepath)) filepath = sub(" \\(.*\\)","",filepath)
        if(!file.exists(filepath)) filepath = paste0("data/",filepath)
        
        # load data and perform initial FFT
        newSpec = readFrameFromFITS(filepath)
        names(newSpec) = tolower(names(newSpec))
        f   = fft(newSpec$flux)
        
        # define attenuation function using double exponential
        att = Vectorize(function(x){1/(exp(exp(input$strength*(x-input$thresh))))})
        
        # calculate weights using attenuation function
        wts   = att(1:length(f))
        
        # rescale according to weights
        f.att = (f * wts)/sum(wts)
        
        # invert FFT, then take real part
        newSpec$f.flt = Re(fft(f.att,inverse=T))
        
        # roughly rescale coefficients using linear regression on flux quantiles
        rescale_coefs = coef(lm(
            quantile(newSpec$flux,probs=c(seq(.15,.85,.05))) ~
                quantile(newSpec$f.flt,probs=c(seq(.15,.85,.05)))
        ))
        newSpec$f.flt = newSpec$f.flt * rescale_coefs[2] + rescale_coefs[1]
        
        # update declared reactive values with new results
        spectrum(newSpec)
        steps(list(f=f,att=att))
        unsaved(TRUE)
    })
    
    # plot periodogram and attenuation curve
    observeEvent(c(input$filepath,input$thresh,input$strength),{ 
        output$freq <- renderPlot({
            ggplot(data.frame(mod=Mod(steps()$f)),aes(x=1:length(mod),y=mod)) +
                geom_line(aes(color="a")) + scale_y_log10(labels = function(x){signif(x,1)}) +
                geom_smooth(aes(color="b"),se=F,method="loess",formula='y~x',span=.1) +
                stat_function(aes(color="c"),fun=function(x){steps()$att(x)+.1},size=1) +
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
    })
    
    # plot original spectrum
    observeEvent(c(input$filepath,input$scale),{
        output$cb58 <- renderPlot({
            switch(input$scale,
                   linear={
                       ggplot(spectrum(),aes(x=loglam,y=flux)) + geom_line() +
                           scale_y_continuous() +
                           labs(x="Log(Wavelength)",y="Flux",
                                title="Original spectrum")
                   },
                   log={
                       ggplot(spectrum(),aes(x=loglam,y=flux)) + geom_line() +
                           scale_y_log10(limits=c(exp(min(log(spectrum()$f.flt),na.rm=T)),max(spectrum()$f.flt))) +
                           labs(x="Log(Wavelength)",y="Log(Flux)",
                                title="Original spectrum")
                   })
            
            
        })
    })
    
    # plot filtered spectrum
    observeEvent(c(input$filepath,input$thresh,input$strength,input$scale),{
        output$cb58flt <- renderPlot({
            switch(input$scale,
                   linear={
                       ggplot(spectrum(),aes(x=loglam,y=f.flt)) + geom_line() +
                           labs(x="Log(Wavelength)",y="Flux",
                                title="Smoothed spectrum")
                   },
                   log={
                       ggplot(spectrum(),aes(x=loglam,y=f.flt)) + geom_line() +
                           scale_y_log10() +
                           labs(x="Log(Wavelength)",y="Log(Flux)",
                                title="Smoothed spectrum")
                   })
        })
    })
    
    # plot noise spectrum
    observeEvent(c(input$filepath,input$thresh,input$strength),{
        output$noise <- renderPlot({
            ggplot(spectrum(),aes(x=loglam,y=flux-f.flt)) + geom_line() +
                labs(x="Log(Wavelength)",y="Flux",
                     title="Noise (original\u2009\u2013\u2009smoothed)")
        })
    })
    
    # save spectrum as template
    observeEvent(input$save,{
        if(unsaved()){
            template = spectrum()
            save(template,file="template.Rdata")
            showNotification("Filtered cB58 spectrum saved as template!",duration=3,type="message")
            unsaved(FALSE)
        }
    })
})
