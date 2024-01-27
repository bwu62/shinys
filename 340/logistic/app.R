library(shiny)

ui = fluidPage(
  titlePanel("Pima logistic"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("model","Which model",c("Null"="null","Glu"="glu","All"="all")),
      sliderInput("threshold","Prediction threshold",0,1,.5,.01),
      verbatimTextOutput("confusionMatrix"),
      verbatimTextOutput("statistics")
    ),
    mainPanel(plotOutput("rocPlot"))
  )
)

server = function(input, output) {
  
  library(MASS)
  library(tidyverse)
  library(pROC)
  library(glue)
  
  if(is.factor(Pima.te$type)) Pima.te$type = ifelse(Pima.te$type=='Yes',1,0)
  model = list(null=glm(type~1, data=Pima.te, family="binomial"),
               glu=glm(type~glu, data=Pima.te, family="binomial"),
               all=glm(type~., data=Pima.te, family="binomial"))
  
  observeEvent(c(input$threshold,input$model),{
    conf.tab = table(
      Pima.te$type,
      factor(predict(model[[input$model]],type="r")>=input$threshold,levels=c(F,T)))
    colnames(conf.tab) = c("Predict 0","Predict 1")
    rownames(conf.tab) = c("No diabetes","Yes diabetes")
    
    output$confusionMatrix = renderPrint({
      conf.tab
    })
    
    output$statistics = renderPrint({
      TP = conf.tab[2,2]
      TN = conf.tab[1,1]
      FP = conf.tab[1,2]
      FN = conf.tab[2,1]
      cat("      Accuracy:",round((TP+TN)/(TP+TN+FP+FN),3))
      cat("\n   Sensitivity:",round(TP/(TP+FN),3))
      cat("\n   Specificity:",round(TN/(TN+FP),3))
      cat("\nPrecision(PPV):",round(TP/(TP+FP),3))
      cat("\n           NPV:", round(TN/(TN+FN),3))
    })
    
    output$rocPlot = renderPlot({
      predictions = predict(model[[input$model]], newdata=Pima.te, type="response")
      ggroc(test_roc <- roc(Pima.te$type~predictions),legacy.axes=T) + 
        labs(title="ROC curve",subtitle=glue("AUC: {round(test_roc$auc,3)}"),
             x="False positive rate",y="True positive rate") + 
        scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
        geom_point(aes(x=conf.tab[1,2]/sum(conf.tab[1,]),y=conf.tab[2,2]/sum(conf.tab[2,])),
                   col="red",size=5) + 
        geom_abline() + coord_equal()
    })
    
  })
  
}

shinyApp(ui = ui, server = server)
