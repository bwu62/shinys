library(shiny)

JScode <-
	"$(function() {
setTimeout(function(){
var vals = [];
for (i = 0; i < 41; i++) {
var val = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,90,100,150,200,250,300,350,400,600,800,1000,2000,5000,10000][i];
val = parseFloat(val.toFixed(8));
vals.push(val);
}
$('#size').data('ionRangeSlider').update({'values':vals})
}, 10)});

var dimension = [0,0];
$(document).on('shiny:connected', function(e) {
dimension[0] = Math.floor(window.innerHeight*0.9);
dimension[1] = Math.floor(window.innerWidth*0.65);
Shiny.onInputChange('dimension', dimension);
});
$(window).resize(function(e) {
dimension[0] = Math.floor(window.innerHeight*0.9);
dimension[1] = Math.floor(window.innerWidth*0.65);
Shiny.onInputChange('dimension', dimension);
});
"

# Define UI for application that draws a histogram
shinyUI(fluidPage(
	
	# Application title
	titlePanel("Central Limit Theorem"),
	
	# Sidebar with a slider input for number of bins 
	sidebarLayout(
		sidebarPanel(
			tags$head(tags$script(HTML(JScode))),
			sliderInput("size",
									"Sample size:",
									min = 1,
									max = 50,
									value = 9),
			radioButtons("dist",
									 "Population distribution:",
									 c("Normal"="normal",
									 	"Bimodal"="bimodal",
									 	"Skewed"="skewed",
									 	"Uniform"="uniform",
									 	"Weibull"="weibull",
									 	"Beta"="beta",
									 	"Sawtooth"="sawtooth"))
		),
		
		# Show a plot of the generated distribution
		mainPanel(
			plotOutput("distPlot",height="auto")
		)
	)
))
