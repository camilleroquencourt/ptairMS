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
ui <- navbarPage("ptairMS",
   
   # Application title
   #titlePanel("ptairMS"),
   
   #tabsetPanel(type = "tabs",
               tabPanel("Read and check data", 
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Create a PtrSet from directory"),
                            shinyFiles::shinyDirButton('dir', 'Select directory', 'Please select a directory',FALSE),
                            shiny::verbatimTextOutput("directory",placeholder = TRUE),
                            shiny::numericInput("mzCalib","mz Reference for calibration",21.022,0,500,1),
                            shiny::actionButton("createPtrSet", "Create PtrSet"),
                            shiny::h3("Load a PtrSet"),
                            shiny::fileInput("ptrSetFile", "Choose RData"),
                            shiny::actionButton("updatePtrSet", "Update PtrSet")
                            ),
                        
                        mainPanel(plotOutput("plotPtrSet")))
                        ),

               tabPanel("Detect peak",
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Detect peak"),
                            shiny::numericInput("PPM","min ppm sepration",130,0,200,10),
                            actionButton("detectPeak", "detect")
                          ),
                          
                          mainPanel(shiny::verbatimTextOutput("peakTable")))),
   
   
               tabPanel("Align samples",
                        # alignement
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::fluidRow(
                              shiny::h3("Alignment between samples"),
                              shiny::selectInput("quantiUnit","Unit:", c("cps", "ncps", "ppb"),selectize=FALSE),
                              shiny::numericInput("fracGroup","Maximum percentage of missing values per features",value = 0.3,min = 0,max = 1,step = 0.1),
                              shiny::numericInput("fracExp",
                                                  "Minimum percentage of samples such that features should be significativly 
                                                  diffrent from baseline:",value = 0.1,min = 0,max = 1,step = 0.1),
                              
                              shiny::column(width = 4, shiny::actionButton('align', 'Align samples')),
                              shiny::column(width = 4,shiny::verbatimTextOutput("nbPeak",placeholder = TRUE)),
                              shiny::column(width = 4,shiny::actionButton('export',"Export")),
                              shiny::column(width = 12,shiny::selectInput("Transformation","Matrix transformation:",c("log2","Centred reduced","None"),
                                                                          selected = FALSE)))
                            
                            ),
                          
                          # Show a plot of the generated distribution
                          mainPanel( plotly::plotlyOutput("table"))
                          
                        )),
   tabPanel("Statistical analysis",
            shiny::fluidRow(
              sidebarPanel(
                shiny::h3("Principal Component Analysis"),
                shiny::numericInput("x","Principal component of x-axis",1,1,4,1),
                shiny::numericInput("y","Principal component of y-axis",2,1,4,1)
              ),
              
              # Show a plot of the generated distribution
              mainPanel(  plotly::plotlyOutput("acp"))
              
            ),
            shiny::fluidRow(
              sidebarPanel(
                shiny::h3("Loading: features order by contribution"),
                shiny::numericInput("maxRow","Number of row to display",10,1,100,1)
              ),
              mainPanel(tableOutput("getLoading"))
            )
            
   )
   #)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rv <- shiny::reactiveValues( ptrset=NULL,path=NULL)
  
  shinyFiles::shinyDirChoose(input, 'dir', roots=c(wd="C:/Users/CR258086/data_copy/voc-transplant/1_raw/ptrms/2019/01 - Janvier/"))
  dir <- reactive(input$dir)
  # path <- reactive({
  #   file.path("C:/Users/CR258086/data_copy",paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  # })

  
  output$directory <- renderPrint({
    if(is.list(input$dir))
      rv$path<-file.path("C:/Users/CR258086/data_copy/voc-transplant/1_raw/ptrms/2019/01 - Janvier/",paste(unlist(input$dir$path[-1]), collapse = .Platform$file.sep))
    rv$path  
    })
  
  shiny::observeEvent(input$createPtrSet, {
   ptrSet<-createPtrSet(rv$path , setName = "ptrSet")
    rv$ptrset<-ptrSet
    })
  
  output$plotPtrSet <- renderPlot({
    if(!is.null(rv$ptrset)) plot(rv$ptrset)
  })
  
  #  shiny::observeEvent(input$close, {
  #    shiny::stopApp(ptrSet)
  #  })
  #  
   shiny::observeEvent(input$detectPeak, {
     ptrSet<-rv$ptrset
     ptrSet<-detectPeak(ptrSet,mz=c(21,29,45,59,69))
     output$peakTable<-renderPrint(ptrSet@peakListAligned)
     rv$ptrset<-ptrSet
   })
  #  
   shiny::observeEvent(input$align, {

     ptrSet<-rv$ptrset
     eSet<- ptairMS::alignSamples(ptrSet)
     
     if(!is.null(eSet)){
       
       output$nbPeak <- renderText({
           X<-Biobase::exprs(eSet)
           paste(nrow(X),"peaks aligned")
       })
       
       output$table<- plotly::renderPlotly({
         
         X<-Biobase::exprs(eSet)
         
         if(input$Transformation == "log2"){
           X<-log2(X)
         } else if (input$Transformation == "Centred reduced"){
           X<-t(apply(X,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
         }
         #ropls::imageF(log2(Biobase::exprs(rv$eset)))
         plotly::plot_ly(z = X,
                         x = colnames(Biobase::exprs(eSet)),
                         y =  paste0("mz-",as.character(rownames(Biobase::exprs(eSet)))), type = "heatmap", colorscale='Jet')
         
       })
     }
  
   })
   
   # nb peak

  
}

#ptrSet <- shiny::runApp(list(ui=ui,server=server))


# Run the application 
shinyApp(ui = ui, server = server)


