#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("ptairMS"),
   
   tabsetPanel(type = "tabs",
               tabPanel("createPtrSet", 
                        fluidRow(shinyFiles::shinyDirButton('dir', 'Directory select', 'Please select a directory',FALSE)),
                        fluidRow(textOutput("caption")),
                        fluidRow(actionButton("createPtrSet", "create")),
                        mainPanel(plotOutput("plotPtrSet")),
                        fluidRow(actionButton("close", "close app"))
                        ),
               tabPanel("detectPeak",
                        fluidRow(actionButton("detectPeak", "detect")),
                        fluidRow(textOutput("peakTable")),
                        fluidRow(actionButton("close2", "close appli"))),
               tabPanel("alignSamples")
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  shinyFiles::shinyDirChoose(input, 'dir', roots=c(wd="C:/Users/CR258086/data_copy/"))
  dir <- reactive(input$dir)
  path <- reactive({
    file.path("C:/Users/CR258086/data_copy",paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  
  output$caption <- renderPrint(path())
  shiny::observeEvent(input$createPtrSet, {
   ptrSet<-createPtrSet(path(),setName = "ptrSet")
   output$plotPtrSet <- renderPlot({
    plot(ptrSet)
   })
   shiny::observeEvent(input$close, {
     shiny::stopApp(ptrSet)
   })
   shiny::observeEvent(input$detectPeak, {
     ptrSet<-detectPeak(ptrSet,mz=59)
     shiny::observeEvent(input$close2, {
       shiny::stopApp(ptrSet)
     })
     output$peakTable<-renderPrint(ptrSet)
   })
   
   
  })
  
}

ptrSet <- shiny::runApp(list(ui=ui,server=server))


# Run the application 
shinyApp(ui = ui, server = server)


