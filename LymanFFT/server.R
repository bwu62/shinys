library(shiny)
library(ggplot2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    # turn off annoying warnings
    options(warn=-1)
    
    # load data and perform initial FFT
    load("./cb58.Rdata")
    f   = fft(cb58$FLUX)
    
    observeEvent(c(input$thresh,input$strength,input$scale),{
        
        # define attenuation function using double exponential
        att = Vectorize(function(x){1/(exp(exp(input$strength*(x-input$thresh))))})
        
        # calculate weights using attenuation function
        wts   = att(1:length(f))
        
        # rescale according to weights
        f.att = (f * wts)/sum(wts)
        
        # invert FFT, then take real part
        cb58$f.flt = Re(fft(f.att,inverse=T))
        
        # roughly rescale coefficients using linear regression on flux quantiles
        rescale_coefs = coef(lm(
            quantile(cb58$FLUX,probs=c(seq(.15,.85,.05))) ~
                quantile(cb58$f.flt,probs=c(seq(.15,.85,.05)))
        ))
        cb58$f.flt = cb58$f.flt * rescale_coefs[2] + rescale_coefs[1]
        
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
                labs(x="Index",y="Relative intensity",title="Fourier transformed spectrum + attenuation curve")
        })
        
        # plot original spectrum
        output$cb58 <- renderPlot({
            ggplot(cb58,aes(x=LOGLAM,y=FLUX)) + 
                geom_line() + switch(input$scale,
                                     linear={scale_y_continuous()},
                                     log={scale_y_log10(
                                         limits=quantile(cb58$FLUX,probs=c(0.01,0.99),na.rm=T)
                                     )})
        })
        
        # plot filtered spectrum
        output$cb58flt <- renderPlot({
            ggplot(cb58,aes(x=LOGLAM,y=f.flt)) + 
                geom_line() + switch(input$scale,
                                     linear={NULL},
                                     log={scale_y_log10()})
        })
        
        
    })
})
