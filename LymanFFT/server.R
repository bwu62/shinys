library(shiny)
library(ggplot2)
library(FITSio)
library(gridExtra)
library(DescTools)

# create server
shinyServer(function(input, output, session) {
    
    # turn off annoying warnings
    options(warn=-1)
    
    # declare reactive values
    steps    = reactiveVal(NA)
    spectrum = reactiveVal(NA)
    unsaved  = reactiveVal(NA)
    initAB   = reactiveVal(TRUE)
    
    # define lm-based rescaling algorithm
    lm.scale = function(source,target){
        breaks = seq(.15,.85,.05)
        return(unname(coef(lm(quantile(target,breaks,na.rm=T)~quantile(source,breaks,na.rm=T)))))
    }
    
    # define loess diff peak extraction algorithm
    loess.peak = function(y,span=0.2,flip=TRUE){
        x=1:length(y)
        return((2*(.5-flip))*(y-predict(loess(y~x,span=span),data.frame(x=x))))
    }
    
    # define low.pass function
    low.pass = function(flux,thresh,strength,debug=FALSE){
        if(debug){
            
            # perform FFT
            f = fft(flux)
            
            # define attenuation function using double exponential
            att = Vectorize(function(x){1/(exp(exp(((1+strength)^2-1)*(x-thresh))))})
            
            # calculate weights using attenuation function
            wts   = att(1:length(flux))
            
            # rescale according to weights, invert, then take real part
            f.flt = Re(fft(f*wts,inverse=T))
            
            # roughly rescale coefficients using linear regression on flux quantiles
            rescale.coefs = lm.scale(f.flt,flux)
            f.flt = f.flt * rescale.coefs[2] + rescale.coefs[1]
            
            # return new values and intermediate steps
            return(list(f.flt=f.flt,steps=list(f=f,att=att)))
            
        } else{
            f.raw = Re(fft(fft(flux)*exp(
                -exp(((1+strength)^2-1)*((1:length(flux))-thresh))
            ),inverse=T))
            ab = lm.scale(f.raw,flux)
            return(f.raw*ab[2]+ab[1])
        }
    }
    
    # run analysis on file
    observeEvent(c(input$filepath,input$thresh,input$strength),{
        
        # find file
        filepath = input$filepath
        if(!file.exists(filepath)) filepath = sub(" \\(.*\\)","",filepath)
        if(!file.exists(filepath)) filepath = paste0("data/",filepath)
        
        # load data
        newSpec = readFrameFromFITS(filepath)
        names(newSpec) = tolower(names(newSpec))
        newSpec = newSpec[,c("loglam","flux")]
        
        # run low pass filter
        res = low.pass(newSpec$flux,input$thresh,input$strength,debug=T)
        
        # update declared reactive values with new results
        newSpec$f.flt = res$f.flt
        spectrum(newSpec)
        steps(res$steps)
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
                labs(x="Frequency index",y="Relative intensity",
                     title="Fourier transformed spectrum + attenuation curve")
        })
    })
    
    # plot original spectrum
    observeEvent(c(input$filepath,input$scale),{
        output$cb58 <- renderPlot({
            switch(input$scale,
                   linear={
                       ggplot(spectrum(),aes(x=loglam,y=flux)) + geom_line() +
                           scale_y_continuous(limits=c(min(spectrum()$f.flt),max(spectrum()$f.flt))) +
                           labs(x="Log(Wavelength)",y="Flux",
                                title="Original spectrum")
                   },
                   log={
                       ggplot(spectrum(),aes(x=loglam,y=flux)) + geom_line() +
                           scale_y_log10(limits=c(exp(min(log(spectrum()$f.flt),na.rm=T)),
                                                  max(spectrum()$f.flt))) +
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
                scale_y_continuous(limits=with(spectrum(),quantile(flux-f.flt,probs=c(.01,.99)))) + 
                labs(x="Log(Wavelength)",y="Flux",
                     title="Noise (original\u2009\u2013\u2009smoothed)")
        })
    })
    
    # save spectrum as template
    observeEvent(input$save,{
        if(unsaved()){
            template = spectrum()[,c("loglam","f.flt")]
            names(template) = c("loglam","flux")
            save(template,file="template.Rdata")
            showNotification("Filtered cB58 spectrum saved as template!",duration=3,type="message")
            unsaved(FALSE)
        }
    })
    
    # plot peak diff curves
    observeEvent(c(input$save,input$filepath3,input$thresh,input$strength,input$span),{
        
        # find file
        filepath3 = input$filepath3
        if(!file.exists(filepath3)) filepath3 = sub(" \\(.*\\)","",filepath3)
        if(!file.exists(filepath3)) filepath3 = paste0("data/",filepath3)
        
        peakSpec = readFrameFromFITS(filepath3)
        names(peakSpec) = tolower(names(peakSpec))
        peakSpec = peakSpec[,c("loglam","flux")]
        peakSpec$f.flt = low.pass(peakSpec$flux,input$thresh,input$strength)
        peakSpec$f.diff = loess.peak(peakSpec$f.flt,span=input$span)
        output$loess = renderPlot({
            ggplot(peakSpec,aes(x=loglam)) +
                geom_line(aes(y=f.flt,color="a")) +
                geom_smooth(aes(y=f.flt,color="b"),method="loess",formula='y~x',
                            span=input$span,size=.5,se=F) +
                scale_color_manual("Legend",
                                   values=c("a"="black","b"="red"),
                                   labels=c("Spectrum","Loess")) +
                theme(legend.justification = c(1, 1), legend.position = c(.98, .98)) +
                labs(x="Log(Wavelength)",y="Flux",
                     title="Smoothed spectrum + loess (for differencing)")
        })
        output$diff = renderPlot({
            ggplot(peakSpec,aes(x=loglam,y=f.diff)) + geom_line() +
                labs(x="Log(Wavelength)",y="Flux difference",
                     title="Trough extraction (loess\u2009\u2013\u2009spectrum)")
        })
        
        if(file.exists("template.Rdata")){
            output$peak.msg = NULL
            load("template.Rdata")
            
            template$f.diff = loess.peak(template$flux,span=input$span)
            conv = data.frame(
                idx = 1:(nrow(peakSpec)-nrow(template)+1),
                raw = convolve(peakSpec$f.diff,template$f.diff,type="filter")
            )
            max_frac = 1/8
            conv$trans = with(conv,pmax(0,abs(raw)-max(raw)*max_frac)^2)
            
            if(startsWith(input$filepath3,"cB58")){
                output$conv = output$trans = output$stats = NULL
            } else{
                output$conv = renderPlot({
                    ggplot(conv,aes(x=idx,y=raw)) +
                        geom_line() +
                        geom_vline(aes(xintercept=which.max(raw)),
                                   color="red",size=1,linetype=2) +
                        geom_ribbon(aes(ymin=-max(raw)*max_frac,max=max(raw)*max_frac),fill="red",alpha=.2) +
                        labs(x="Wavelength index",y="Convolution",
                             title="Convolution with trough-extracted template")
                })
                
                output$trans = renderPlot({
                    ggplot(conv,aes(x=idx,y=trans)) +
                        geom_line() +
                        # geom_vline(aes(xintercept=which.max(conv)),
                        #            color="red",size=1,linetype=2) +
                        labs(x="Wavelength index",y="Score",
                             title="Transformed convolution (to emphasize Gini)")
                })
                
                output$stats = renderTable(
                    data.frame(
                        vars = c("Index of max:","Raw Gini index:","Transformed Gini:"),
                        vals = with(conv,as.character(
                            c(which.max(trans),round(Gini(abs(raw)),3),round(Gini(trans),3))
                        ))
                    ),
                    spacing="xs",align="l",colnames=FALSE
                )
            }
        } else{
            output$peak.msg = renderText(paste0("Cannot match peaks without template; ",
                                                "save smoothed template first!"))
            output$conv = NULL
        }
    })
    
    observeEvent(input$offset,{
        initAB(TRUE)
    })
    
    # # plot offset
    # # input = list(filepath2="spec-5328-55982-0218.fits",thresh=100,
    # #              strength=.03,stretch=2.9,shift=-.5,offset=530)
    # observeEvent(c(input$save,input$filepath2,input$offset,input$stretch,input$shift),{
    #     if(file.exists("template.Rdata")){
    #         output$match.msg = NULL
    #         load("template.Rdata")
    #         matchSpec = readFrameFromFITS(paste0("data/",input$filepath2))[,c("loglam","flux")]
    #         names(matchSpec) = tolower(names(matchSpec))
    # 
    #         if(nrow(matchSpec) < nrow(template)){
    #             output$offset = NULL
    #             output$match.msg = renderText("Input spectrum shorter than template!")
    #         } else{
    #             matchSpec$f.flt = low.pass(matchSpec$flux,input$thresh,input$strength)
    #             N = nrow(matchSpec)-nrow(template)
    #             updateSliderInput(session,"offset",max=N)
    #             matchSpec$tempflux = c(rep(NA,input$offset),template$flux,rep(NA,N-input$offset))
    #             
    #             matchSpec.complete = matchSpec[complete.cases(matchSpec),]
    #             
    #             if(initAB()){
    #                 ab = with(matchSpec.complete,lm.scale(tempflux,f.flt))
    #                 updateSliderInput(session,"stretch",value=ab[2])
    #                 updateSliderInput(session,"shift",value=ab[1])
    #                 initAB(FALSE)
    #             }
    #             
    #             output$offset = renderPlot({
    #                 ggplot(matchSpec,aes(x=loglam)) +
    #                     geom_line(aes(y=f.flt,color="a")) +
    #                     geom_line(aes(y=tempflux*input$stretch+input$shift,color="b")) + 
    #                     scale_color_manual("Legend",
    #                                        values=c("a"="red","b"="black"),
    #                                        labels=c("Spectrum","Template")) +
    #                     theme(legend.justification = c(1, 1), legend.position = c(.98, .98)) 
    #             })
    #             
    #             output$quants = renderPlot({
    #                 breaks = seq(.15,.85,.05)
    #                 grid.arrange(
    #                     # ggplot(matchSpec.complete) + 
    #                     #     stat_ecdf(aes(f.flt,color="a")) + 
    #                     # stat_ecdf(aes(tempflux*input$stretch+input$shift,color="b"))  + 
    #                     # scale_color_manual("Legend",
    #                     #                    values=c("a"="red","b"="black"),
    #                     #                    labels=c("Spectrum","Template")) +
    #                     # theme(legend.justification = c(1, 0), legend.position = c(.98, .02)) , 
    #                     ggplot(data.frame(lapply(matchSpec.complete,
    #                                              quantile,probs=breaks,na.rm=T))) + 
    #                         geom_point(aes(x=breaks,y=f.flt)) + 
    #                         geom_point(aes(x=breaks,y=tempflux*input$stretch+input$shift)) ,
    #                     ggplot(matchSpec.complete,aes(x=loglam,y=f.flt-tempflux)) + geom_point() , 
    #                     ncol=2
    #                 )
    #             })
    #         }
    #     } else{
    #         output$match.msg = renderText("Save smoothed template first!")
    #         output$offset = output$ecdf = NULL
    #     }
    # })
})
