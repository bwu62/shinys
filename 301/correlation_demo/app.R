library(shiny)

ui <- fluidPage(
  
  titlePanel("Correlation demo"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("r",
                  "Correlation coefficient",
                  min = -1,
                  max = 1,
                  value = 0,
                  step = .01)
    ),
    
    mainPanel(
      plotOutput("points")
    )
  )
)

server <- function(input, output) {
  
  data = reactive(data.frame(mvtnorm::rmvnorm(1000,mean = c(0,0),sigma = matrix(c(1,input$r,input$r,1),2))))
  
  output$points <- renderPlot({
    par(las=1,pty="s",cex=.7)
    with(data(),plot(X1,X2,xlim=c(-4,4),ylim=c(-4,4),pch=20))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
