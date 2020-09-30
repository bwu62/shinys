library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    titlePanel("Lyman FFT low-pass smoothing"),
    
    sidebarLayout(
        sidebarPanel(
            tags$head(tags$style("#freq{height:25vh !important;}
                                 #cb58{height:25vh !important;margin-top:2vh !important;}
                                 #cb58flt{height:25vh !important;margin-top:2vh !important;}
                                 #noise{height:25vh !important;margin-top:2vh !important;}")),
            selectInput("filepath","Choose input spectrum:",
                        c("cB58_Lyman_break.fit",list.files("data/"))),
            sliderInput("thresh",
                        "Filter threshold:",
                        min   = 0,
                        max   = 500,
                        value = 100,
                        step  = 10),
            
            sliderInput("strength",
                        "Attenuation strength:",
                        min   = 0.000,
                        max   = 0.05,
                        value = 0.02,
                        step  = 0.001),
            radioButtons("scale",
                         "Spectra flux scale:",
                         c("linear","log"))
        ),
        mainPanel(
            fluidPage(
                plotOutput("freq"),
                plotOutput("cb58"),
                plotOutput("cb58flt"),
                plotOutput("noise")
            )
        )
    )
))
