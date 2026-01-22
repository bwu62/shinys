library(shiny)

ui = fluidPage(
  
  titlePanel("Tint laws"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("type","Type of vehicle",c("Car"="car","SUV/Truck"="mpv")),
      sliderInput("FRONT","Front side window",0,100,100,5),
      sliderInput("BACK","Back side window",0,100,100,5),
      sliderInput("REAR","Rear window",0,100,100,5),
      radioButtons("WIND","Windshield",c(
        "None"="none","4-6 inches OR AS-1 line"="as1_46","70% or higher VLT"="vlt70")),
      
    ),
    
    mainPanel(
      plotOutput("legality")
    )
  )
)
# 
server = function(input, output) {
  
  library(ggplot2) ; library(dplyr) ; library(scales) ; library(glue)
  load("tint_data.Rdata")
  
  output$legality = renderPlot({
    
    tints.us = spdf_fortified[[input$type]] %>% 
      mutate(legal = factor((input$FRONT/100>=front_side) & (input$BACK/100>=back_side) & 
               (input$REAR/100>=rear) & 
               (input$WIND=="none" | (input$WIND=="as1_46" & windshield!="No tint allowed") | 
                  (input$WIND=="vlt70" & windshield=="No more than 70%")),
               levels=c(TRUE,FALSE),labels=c("Legal","Illegal")))
    
    ggplot(tints.us) + 
      geom_polygon(aes(x=long,y=lat,group=state,fill=legal), color="gray50") +
      geom_text(data=centers[[input$type]],aes(x,y,label=id),color="white",size=4) +
      scale_fill_manual(values=c("#00bfc4","#f8766d"),labels=c("Legal","Illegal"),drop=F) +
      ggtitle(glue("Your vehicle is legal in {sum(2-as.numeric(tints.us$legal))/7}/51 or ",
      "{round(100*mean(2-as.numeric(tints.us$legal)),0)}% of states(+DC)")) + 
      theme_void() + coord_map()
    
  })
}


shinyApp(ui=ui,server=server)
