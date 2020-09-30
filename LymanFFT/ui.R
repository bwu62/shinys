library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    titlePanel("Lyman FFT low-pass smoothing"),

    sidebarLayout(
        sidebarPanel(
            tags$head(tags$style("#freq{height:29vh !important;}
                                 #cb58{height:29vh !important;margin-top:2vh !important;}
                                 #cb58flt{height:29vh !important;margin-top:2vh !important}")),
            sliderInput("thresh",
                        "Filter threshold:",
                        min = 0,
                        max = 500,
                        value = 100,
                        step=10),
            sliderInput("strength",
                        "Attenuation strength:",
                        min = 0.000,
                        max = 0.05,
                        value = 0.03,
                        step=0.001),
            radioButtons("scale",
                         "Spectrum flux scale:",
                         c("linear","log"))
        ),
        mainPanel(
            fluidPage(
                plotOutput("freq",height="30%"),
                plotOutput("cb58",height="30%"),
                plotOutput("cb58flt",height="30%")
            )
        )
    )
))
