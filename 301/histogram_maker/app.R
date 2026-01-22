library(shiny)

ui <- fluidPage(
  titlePanel("Histogram maker"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        inputId = "text",
        label = "Enter data (separated by spaces, commas, or new lines)"),
      textInput(
        inputId = "bins",
        label = paste0("Enter class interval boundaries (separated by commas). ",
        "If empty, or if does not completely cover data, bins will be automatically chosen.")),
      radioButtons(
        inputId = "freq",
        label = "What type of scale?",
        c("Density"=FALSE,"Frequency"=TRUE)
      ),
      radioButtons(
        inputId = "right",
        label = "Include which endpoint of each class interval?",
        c("Left"=FALSE,"Right"=TRUE)
      ),
      textInput(
        inputId = "title",
        label = "Title:", value = "Histogram of ___"
      ),
      textInput(
        inputId = "xlab",
        label = "Horizontal label:", value = "Name of variable"
      )
    ),
    mainPanel(
      plotOutput("hist")
    )
  )
)

server <- function(input, output) {
  
  data = reactive(
    readr::parse_number(
      stringr::str_subset(
        stringr::str_split(input$text,"\\n|,| ")[[1]],
        ".+"
      )
    )
  )
  
  bins = reactive({
    temp = try(
      readr::parse_number(
        stringr::str_subset(
          stringr::str_split(input$bins,",| ")[[1]],
          ".+"
        )
      )
    )
    if(length(temp)<2||!is.numeric(temp)||max(temp)<max(data())||min(temp)>min(data())){
      "Sturges"
    } else{
      temp
    }
  })
  
  output$hist <- renderPlot({
    try(
      if(input$freq){
        hist(data(),breaks=bins(),right=input$right,main=input$title,xlab=input$xlab)
      } else{
        hist(data(),breaks=bins(),freq=FALSE,right=input$right,main=input$title,xlab=input$xlab)
      }
    )
  })
}

shinyApp(ui = ui, server = server)
