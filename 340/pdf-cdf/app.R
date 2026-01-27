library(shiny)
library(gridExtra)
library(tidyverse)

ui <- fluidPage(

    titlePanel("PDF, CDF, inverse CDF for Normal"),

    sidebarLayout(
        sidebarPanel(
          sliderInput("mu", label="Mean",min=-3,max=3,value=0,step=0.1),
          sliderInput("sigma2", label="Variance",min=.01,max=10,value=1,step=.1),
          sliderInput("x", label="x",min=-3,max=3,value=0,step=0.1)
        ),

        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      grid.arrange(
        ggplot()+stat_function(fun=~{dnorm(.x,input$mu,sqrt(input$sigma2))},n=1001)+
          scale_x_continuous(expand=c(0,0),limits=c(-3,3)) + 
          scale_y_continuous(expand=c(.005,.005)) + 
          geom_point(aes(x=input$x,y=dnorm(input$x,input$mu,sqrt(input$sigma2))),color="red",size=4) + 
          stat_function(fun=~{dnorm(.x,input$mu,sqrt(input$sigma2))},n=1001,geom="area",xlim=c(-100,input$x),fill="red",alpha=.5)+
          labs(title="Probability density function (PDF)",x="x",y="Density"),
        ggplot()+geom_function(fun=~{pnorm(.x,input$mu,sqrt(input$sigma2))},n=1001)+
          geom_point(aes(x=input$x,y=pnorm(input$x,input$mu,sqrt(input$sigma2))),color="red",size=4) +
          scale_x_continuous(expand=c(.005,.005),limits=c(-3,3)) + scale_y_continuous(expand=c(.005,0)) + labs(title="Cumulative distribution function (CDF)",x="x",y="Probability"),
        ggplot()+geom_function(fun=~{qnorm(.x,input$mu,sqrt(input$sigma2))},n=1001)+
          geom_point(aes(x=pnorm(input$x,input$mu,sqrt(input$sigma2)),y=input$x),color="red",size=4) +
          scale_x_continuous(expand=c(.005,.005),limits=c(0,1)) + scale_y_continuous(expand=c(.005,.005),limits=c(-3,3)) + labs(title="Quantile function (inverse CDF)",x="Probability",y="x"),
        layout_matrix=matrix(c(1,1,1,3,3,2,2,2,3,3),ncol=5,byrow=T))
    }#,height=600
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
