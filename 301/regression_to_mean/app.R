library(shiny)
library(mvtnorm)
library(tidyverse)

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      label input[value='y']+span{color:blue;font-weight:bold;}
      label input[value='x']+span{color:firebrick;font-weight:bold;}
      div.radio:not(:first-child) label input[type='radio'][name='rtm']+span{color:darkgreen;font-weight:bold;}
  "))),
  
  titlePanel("Regression to mean demo"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("r",
                  "Correlation coefficient",
                  min = -1,
                  max = 1,
                  value = 0.5,
                  step = .01),
      checkboxGroupInput("lines",
                         "Lines to show:",
                         choices = c(
                           "SD line" = "s",
                           "Y on X regression line" = "y",
                           "X on Y regression line" = "x"
                         ),
                         selected = "y"),
      radioButtons("rtm",
                   "Regression to mean effect:",
                   choices = c(
                     "Hide" = 0,
                     "Y on X (vertical slices)" = 1,
                     "X on Y (horizontal slices)" = -1
                   ),
                   selected = 0),
      conditionalPanel(
        condition = "input.rtm != 0",
        sliderInput("slice",
                    "Where to slice?",
                    min = -4,
                    max = 4,
                    value = 0,
                    step = 0.1)
      )
    ),
    
    mainPanel(
      plotOutput("points")
    )
  )
)

server <- function(input, output) {
  
  data = reactive(setNames(data.frame(rmvnorm(2000,mean = c(0,0),sigma = matrix(c(1,input$r,input$r,1),2))),c("x","y")))
  
  output$points <- renderPlot({
    
    slice_width = 0.2
    data = data()
    mx = mean(data$x) ; my = mean(data$y)
    sx = sd(data$x)   ; sy = sd(data$y)
    r = cor(data$x,data$y)
    
    ggplot(data,aes(x,y)) + geom_point(size=.6,alpha=.6) + theme_bw() + coord_equal() +
      geom_abline(slope=sign(r),color=ifelse("s" %in% input$lines,"black",NA)) + 
      geom_abline(slope=r,color=ifelse("y" %in% input$lines,"blue",NA),linewidth=1) +
      geom_abline(slope=1/r,color=ifelse("x" %in% input$lines,"firebrick",NA),linewidth=1) +
      scale_x_continuous(breaks=-4:4,limits=c(-4,4),expand=0,minor_breaks=NULL) +
      scale_y_continuous(breaks=-4:4,limits=c(-4,4),expand=0,minor_breaks=NULL) + 
      geom_vline(xintercept=input$slice+slice_width*c(-1,1),color=ifelse(1==input$rtm,"darkgreen",NA)) + 
      geom_point(data = data %>% filter(abs(x-input$slice)<slice_width),aes(x,y),
                 color=ifelse(1==input$rtm,"darkgreen",NA),size=.7) + 
      geom_hline(yintercept = r*input$slice,
                 color=ifelse(1==input$rtm,"darkgreen",NA),linetype="dotted",linewidth=1) + 
      geom_hline(yintercept=input$slice+slice_width*c(-1,1),color=ifelse(-1==input$rtm,"darkgreen",NA)) + 
      geom_point(data = data %>% filter(abs(y-input$slice)<slice_width),aes(x,y),
                 color=ifelse(-1==input$rtm,"darkgreen",NA),size=.7) + 
      geom_vline(xintercept = 1/r*input$slice,
                 color=ifelse(-1==input$rtm,"darkgreen",NA),linetype="dotted",linewidth=1)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
