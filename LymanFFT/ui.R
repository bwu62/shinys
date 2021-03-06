library(shiny)

shinyUI(
    navbarPage(
        "Lyman-break galaxies",
        tabPanel(
            "Smoothing",
            sidebarLayout(
                sidebarPanel(
                    tags$head(tags$style("#freq{height:23vh !important;}
                                          #cb58 {height:23vh !important;margin-top:0vh !important;}
                                          #cb58flt{height:23vh !important;margin-top:0vh !important;}
                                          #noise{height:23vh !important;margin-top:0vh !important;}")),
                    selectInput("filepath","Choose input spectrum:",
                                c("cB58_Lyman_break.fit (template)",
                                  list.files("data/",pattern=paste0("shapley.*|spec-(1353|5328|6064|5972|2557|2261|5793|",
                                                                    "1568|6008|1437|0390|2124|5324|530.|6184|7258|3808|6439|7257|5454).*")))),
                    sliderInput("thresh",
                                "Filter threshold:",
                                min   = 0,
                                max   = 500,
                                value = 200,
                                step  = 10),
                    
                    sliderInput("strength",
                                "Attenuation strength:",
                                min   = 0.000,
                                max   = 0.05,
                                value = 0.03,
                                step  = 0.001),
                    radioButtons("scale",
                                 "Spectra flux scale:",
                                 c("linear","log")),
                    conditionalPanel(
                        "input.filepath == 'cB58_Lyman_break.fit (template)'",
                        actionButton("save","Save template")
                    ),
                    textOutput("saved")
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
        ),
        tabPanel(
            "Peak matching",
            sidebarLayout(
                sidebarPanel(
                    tags$head(tags$style("#loess{height:23vh !important;}
                                          #diff{height:23vh !important;}
                                          #conv{height:23vh !important;}
                                          #trans{height:23vh !important;}")),
                    selectInput("filepath2","Choose spectrum for peak extraction",
                                c("cB58_Lyman_break.fit (template)",
                                  list.files("data/",pattern=paste0("shapley.*|spec-(1353|5328|6064|5972|2557|2261|5793|",
                                                                    "1568|6008|1437|0390|2124|5324|530.|6184|7258|3808|6439|7257|5454).*")))
                    ),
                    sliderInput("span",
                                "Set loess span",
                                min   = 0.01,
                                max   = 0.5,
                                value = 0.2,
                                step  = 0.01),
                    # span(textOutput("max"),style="font-size:16px;"),
                    # span(textOutput("gini"),style="font-size:16px;")
                    conditionalPanel("output$stats",tableOutput("stats"))
                ),
                mainPanel(
                    plotOutput("loess"),
                    plotOutput("diff"),
                    span(textOutput("peak.msg"),style="font-size:18px;font-weight:bold;color:red;"),
                    conditionalPanel("output$conv",plotOutput("conv")),
                    conditionalPanel("output$trans",plotOutput("trans"))
                )
            )
        ),
        tabPanel(
            "Scale matching",
            sidebarLayout(
                sidebarPanel(
                    tags$head(tags$style("#conv2{height:28vh !important;}
                                          #offset{height:28vh !important;}
                                          #quants{height:30vh !important;}")),
                    span("Showing match for file:",style="font-size:18px;font-weight:bold;"),
                    span(textOutput("file"),style="font-size:14px"),
                    HTML("<br>"),
                    sliderInput("offset",
                                "Template offset:",
                                min   = 0,
                                max   = NA,
                                value = 0,
                                step  = 1),
                    sliderInput("stretch",
                                "Template stretch:",
                                min   = 0.0001,
                                max   = 1000,
                                value = 1,
                                step  = 0.001),
                    sliderInput("shift",
                                "Template shift:",
                                min   = -1000,
                                max   = 1000,
                                value = 0,
                                step  = 0.001),
                    conditionalPanel("output$stats2",tableOutput("stats2"))
                ),
                mainPanel(
                    conditionalPanel("output$match.msg",textOutput("match.msg")),
                    conditionalPanel("output$conv2",plotOutput("conv2")),
                    conditionalPanel("output$offset",plotOutput("offset")),
                    conditionalPanel("output$quants",plotOutput("quants"))
                )
            )
        )
    )
)