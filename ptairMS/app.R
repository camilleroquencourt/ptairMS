#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ptairMS)
library(dplyr)

score_plotly <- function(ropls.model,label.c = "sampleNames",color.c = "",info.vc = "all",
                         title.c = "", components.vi = c(1, 2),palette.c = "Set1",ellipse.l = TRUE,
                         size.ls = list(axis_lab.i = 16,
                                        axis_text.i = 14,
                                        point.i = 3,
                                        label.i = 5,
                                        title.i = 20,
                                        legend_title.i = 15,
                                        legend_text.i = 15),
                         figure.c = c("interactive",
                                      "interactive_plotly",
                                      "my_scoreplot.pdf",
                                      "my_scoreplot.html")[2]) 
{
  
  
  # checking arguments and preparing data
  
  ## dataset [ExpressionSet]
  
  eset <- ropls::getEset(ropls.model)
  
  if (any(dim(eset) < 1))
    stop("ExpressionSet object could not be extracted from the 'ropls.model'")
  
  pdata.df <- Biobase::pData(eset)
  #pdata.df$date<-as.factor(pdata.df$date)
  data.df <- cbind.data.frame(sampleNames = Biobase::sampleNames(eset),
                              pdata.df,
                              stringsAsFactors = FALSE)
  
  ## plotly info
  
  if (length(info.vc) == 1 && info.vc == "all") {
    text.df <- data.df
  } else {
    info_in_metadata.vl <- info.vc %in% colnames(data.df)
    if (any(!info_in_metadata.vl)) {
      if (all(!info_in_metadata.vl)) {
        stop("None of the selected info.vc names was found in the sampleMetadata:\n",
             paste(info.vc, collapse = ", "))
      } else {
        warning("The following columns were not found in the sampleMetadata:\n",
                paste(info.vc[!info_in_metadata.vl], collapse = ", "))
      }
    }
    text.df <- data.df[, info.vc[info_in_metadata.vl], drop = FALSE]
  }
  
  text.vc <- apply(text.df, 1,
                   function(row.vc)
                     paste(paste0(colnames(text.df), " = ", row.vc), collapse = "\n"))
  
  ## score vectors
  
  score.mn <- ropls::getScoreMN(ropls.model)
  data.df <- cbind.data.frame(data.df,
                              text = text.vc,
                              score.mn)
  
  ## labels
  
  stopifnot(length(label.c) == 1)
  
  if (label.c != "")
    stopifnot(label.c %in% colnames(data.df))
  
  ## colors
  
  stopifnot(length(color.c) == 1)
  
  if (color.c != "") {
    stopifnot(color.c %in% colnames(data.df))
    if (is.factor(data.df[, color.c]) || is.character(data.df[, color.c])) {
      color_type.c <- "qualitative"
    } else
      color_type.c <- "quantitative"
  }
  
  ## components
  
  stopifnot(length(components.vi) == 2)
  stopifnot(max(components.vi) <= ncol(score.mn))
  
  ## sizes
  
  size_default.vi <- c(axis_lab.i = 16,
                       axis_text.i = 14,
                       point.i = 3,
                       label.i = 5,
                       title.i = 20,                      
                       legend_title.i = 15,
                       legend_text.i = 15)
  
  for (size.c in names(size_default.vi)) {
    if (!(size.c %in% names(size.ls)))
      size.ls[[size.c]] <- size_default.vi[size.c]
  }
  
  ## filename extention
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c), ".", fixed = TRUE)), 1)
  
  
  # starting the plot [ggplot]
  
  p <- eval(parse(text = paste0("ggplot2::ggplot(data.df, ggplot2::aes(x = ",
                                paste0('p', components.vi[1]),
                                ", y = ",
                                paste0('p', components.vi[2]),
                                ifelse(color.c != "",
                                       paste0(", color = ", color.c),
                                       ""),
                                ifelse(label.c != "",
                                       paste0(", label = ", label.c),
                                       ""),
                                ", text = text))")))
  
  ## text/points [geom_text/geom_point]
  
  if (label.c != "") {
    p <- p + ggplot2::geom_text(size = size.ls[["label.i"]], ggplot2::aes(fontface = "bold"))
  } else {
    p <- p + ggplot2::geom_point(size = size.ls[["point.i"]])
  }
  
  ## horizontal and vertical lines [geom_hline, geom_vline]
  
  p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0))
  
  ## ellipses [stat_ellipse]
  
  p <- p + eval(parse(text = paste0("ggplot2::stat_ellipse(ggplot2::aes(x = ",
                                    paste0('p', components.vi[1]),
                                    ", y = ",
                                    paste0('p', components.vi[2]),
                                    ", group = 1), type = 'norm')")))
  
  if (ellipse.l && color.c != "" && color_type.c == "qualitative")
    p <- p + eval(parse(text = paste0("ggplot2::stat_ellipse(ggplot2::aes(x = ",
                                      paste0('p', components.vi[1]),
                                      ", y = ",
                                      paste0('p', components.vi[2]),
                                      ", group = ",
                                      ifelse(color.c != "", color.c, 1),
                                      "), type = 'norm')")))
  
  # title and axis labels [labs]
  
  p <- p + ggplot2::labs(title = title.c,
                         x = paste0("t", components.vi[1],
                                    " (",
                                    round(ropls.model@modelDF[components.vi[1], "R2X"] * 100),
                                    "%)"),
                         y = paste0("t", components.vi[2],
                                    " (",
                                    round(ropls.model@modelDF[components.vi[2], "R2X"] * 100),
                                    "%)"))
  
  # theme [them_bw, theme]
  
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = size.ls[["title.i"]], face = "bold"),
                   axis.title.x = ggplot2::element_text(size = size.ls[["axis_lab.i"]], face = "bold"),
                   axis.title.y = ggplot2::element_text(size = size.ls[["axis_lab.i"]], face = "bold"),
                   axis.text = ggplot2::element_text(size = size.ls[["axis_text.i"]]),
                   legend.title = ggplot2::element_text(face = "bold", size = size.ls[["legend_title.i"]]),
                   legend.text = ggplot2::element_text(face = "bold", size = size.ls[["legend_text.i"]]))
  
  # palette [scale_colour_brewer, scale_colour_gradientn]
  
  if (color.c != "") {
    if (color_type.c == "qualitative") {
      if (palette.c != "")
        p <- p + ggplot2::scale_colour_brewer(palette = palette.c)
    } else
      p <- p + ggplot2::scale_colour_gradientn(colours = rev(rainbow(100, end = 4/6)))
  }
  
  # display/saving [plotly::ggplotly, plotly::layout, htmlwidgets::saveWidget, plotly::as_widget]
  
  if (figure.c == "interactive_plotly" || filename_ext.c == "html") {
    
    p <- plotly::ggplotly(p, tooltip = "text")
    
    p <- plotly::layout(p, hoverlabel = list(font = list(size = 20)))
    
    if (filename_ext.c == "html") {
      
      htmlwidgets::saveWidget(plotly::as_widget(p), figure.c)
      
      return(p)
      
      
    } else {
      
      #print(p)
      
      return(p)
      
    }
    
  } else {
    
    if (filename_ext.c == "pdf")
      grDevices::pdf(figure.c)
    
    #print(p)
    
    if (filename_ext.c == "pdf")
      grDevices::dev.off()
    
    return(p)
    
  }
  
}               


# Define UI for application that draws a histogram
ui <- navbarPage("ptairMS",
   
   # Application title
   #titlePanel("ptairMS"),
   
   #tabsetPanel(type = "tabs",
               tabPanel("Read and check data", 
                        
                        shiny::fluidPage(
                          shiny::h3("Create/Load a PtrSet from directory"),
                          shiny::selectInput("createOrLoad", "Choose action", c("loadPtrSet","createPtrSet"),
                                      multiple=FALSE),
                          uiOutput("param")
                        ),
                        
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Your PtrSet:"),
                            shiny::verbatimTextOutput("showPtrSet",placeholder = FALSE),
                            uiOutput("update")
                            ),
                        
                        mainPanel( 
                                   plotOutput("plotPtrSet"))),
                        
                        
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Calibration"),
                            uiOutput("paramCalib")
                          ),
                          
                          mainPanel( 
                            plotOutput("plotCalib"))),
                        
                        shiny::sidebarLayout(
                          sidebarPanel(
                           shiny::h3("Expirations phases"),
                            uiOutput("paramTimeLimit")
                          ),
                          mainPanel(                          
                            shiny::splitLayout(
                              cellWidths = c("30%","70%"),
                            DT::dataTableOutput("tableTimeLimit"),
                            plotly::plotlyOutput("plotTimeLimit")
                            
                          ))

                          )
                           

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
                        
                        
                        # alignement
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::fluidPage(
                              shiny::h3("Alignment between samples"),
                              shiny::selectInput("quantiUnit","Unit:", c("ppb", "ncps", "cps"),selectize=FALSE),
                              shiny::numericInput("fracGroup","Maximum percentage of missing values per features",value = 0.3,min = 0,max = 1,step = 0.1),
                              shiny::numericInput("fracExp",
                                                  "Minimum percentage of samples such that features should be significativly 
                                                  diffrent from baseline:",value = 0.1,min = 0,max = 1,step = 0.1),
                              shiny::fluidRow(
                                shiny::column(width = 6, checkboxInput("filterExp", label = "variable with oscillation", value = TRUE))
                              ),
                              shiny::fluidRow(
                              shiny::column(width = 4, shiny::actionButton('align', 'Align samples')),
                              shiny::column(width = 4,shiny::verbatimTextOutput("nbPeak",placeholder = TRUE)),
                              shiny::column(width = 4,shiny::actionButton('export',"Export"))),
                              shiny::selectInput("Transformation","Matrix transformation:",c("log2","Centred reduced","None"),
                                                                         selected = FALSE),
                              shiny::fluidRow(
                              shiny::column(width = 6, shiny::actionButton("impute", label = "Impute")),
                              shiny::column(width = 6,shiny::actionButton("outlier","Delete outlier")))
                              
                              )),
                          
                          # Show a plot of the generated distribution
                          mainPanel( plotly::plotlyOutput("tableAlign"))),

                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Modeled and raw data"),
                            shiny::uiOutput("secondSelectionModeledData")
                          ),
                          
                          # Show a plot of the generated distribution
                          mainPanel(plotly::plotlyOutput("dataModeled"))
                          
                        ),
                        shiny::sidebarLayout(
                          sidebarPanel(
                            shiny::h3("Longitudinal vizualisation"),
                            shiny::uiOutput("secondSelectionLongitudinal")
                          ),
                          
                          # Show a plot of the generated distribution
                          mainPanel(plotly::plotlyOutput("dataLongitudinal"))
                          
                        )
                        
                        ),
   tabPanel("Statistical analysis",
            shiny::fluidRow(
              sidebarPanel(
                shiny::h3("Principal Component Analysis"),
                shiny::numericInput("x","Principal component of x-axis",1,1,4,1),
                shiny::numericInput("y","Principal component of y-axis",2,1,4,1),
                shiny::uiOutput("pca")
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
  
  paramRv<-shiny::reactiveValues(path=NULL)
  rv <- shiny::reactiveValues( ptrset=NULL,dataTimeLimit=NULL,bgPoints=NULL)
  
  ## user interface 
  output$param<-renderUI({
    if(input$createOrLoad == "createPtrSet"){
      shiny::fluidPage(
        shiny::fluidRow(
          shiny::column(3,shiny::actionButton("selectDir","Select a directory")),
          shiny::column(6,shiny::verbatimTextOutput(outputId = "showDirectorySelect",
                                                    placeholder = TRUE)),
          shiny::column(3,shiny::textInput(inputId = "ptrName",label = "set name"))),
        
        shiny::fluidRow(shiny::h4("Calibration masses"),
        shiny::column(4,shiny::numericInput(inputId = "mz1",label = "mz1",value = 21.022)),
        shiny::column(4,shiny::numericInput(inputId = "mz2",label = "mz2",value = 29.013424)),
        shiny::column(4,shiny::numericInput(inputId = "mz3",label = "mz3",value = 60.0525)),
        shiny::column(4,shiny::numericInput(inputId = "mz4",label = "mz4",value = 203.943)),
        shiny::column(4,shiny::numericInput(inputId = "mz5",label = "mz5",value = 330.8495))),
        
        shiny::fluidRow(
          shiny::column(6,shiny::checkboxGroupInput(inputId = "BoolExp",
                                                    label = "Detect expirations phases",
                                                    choiceNames = c("TRUE"), choiceValues = c(TRUE))),
          shiny::column(6,shiny::numericInput(inputId = "mzBreathTracer",label = "mzBreathTracer",value = 44.997))
         ),
        #saveDir
        shiny::fluidRow(
          shiny::p(class = 'text-center',
                   shiny::actionButton(inputId = "createPtrSet",label = "create a PtrSet")
         ))
      )
    } else {
      shiny::fluidRow( 
        shiny::p(class = 'text-center',
          shiny::actionButton("selectFile","Choose a .Rdata ptrSet object"))
        )
    }
  })
  
  shiny::observeEvent(input$selectDir, {
    pth <- choose.dir() #tk_choose.dir
    paramRv$path<-pth
  })
  
  output$showDirectorySelect<-renderPrint({
    if(!is.null(paramRv$path)) print(paramRv$path)
  })
  
  shiny::observeEvent(input$createPtrSet, {
    if(!is.null(input$BoolExp)) fracMaxTIC<-0.7 else fracMaxTIC<-0
    ptrSet<-createPtrSet(dir = paramRv$path , setName = input$ptrName,
                         mzCalibRef = c(input$mz1,input$mz2,input$mz3,input$mz4,input$mz5),
                         fracMaxTIC = fracMaxTIC,
                         mzBreathTracer = input$mzBreathTracer, saveDir = NULL)
    rv$ptrset<-ptrSet
    })
  
  shiny::observeEvent(input$selectFile, {
    pth <- file.choose()
    load_obj <- function(f){
       env <- new.env()
       nm <- load(f, env)[1]
       env[[nm]]
      }
    rv$ptrset<-load_obj(pth)
  })
  
  
  output$showPtrSet<-renderPrint(
    if(!is.null(rv$ptrset)) print(rv$ptrset)
    )
    
  output$plotPtrSet <- renderPlot({
    if(!is.null(rv$ptrset)) plot(rv$ptrset)
  })
  
  output$update<-renderUI({
    if(!is.null(rv$ptrset)) shiny::actionButton(inputId = "updatePtrSet",label = "Update ptrSet")
  })
  
  
  shiny::observeEvent(input$updatePtrSet,{
    if(!is.null(rv$ptrset)){
      rv$ptrset<-updatePtrSet(rv$ptrset)
    }
  })
  

  
  output$paramTimeLimit<-renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(

        shiny::fluidRow(
          shiny::selectInput("fileNameTimeLimit","File name:",as.list(names(rv$ptrset@TIC)),selectize=FALSE)
        ),
        shiny::fluidRow(shiny::h5("parameters")),
        shiny::fluidRow(
            shiny::column(width = 4,shiny::numericInput("fracMaxTIC", "FracMaxTIC", rv$ptrset@parameter$timeLimit$fracMaxTIC,0, 1, 0.05)),
            shiny::column(width = 4,shiny::numericInput("fracMaxTICBg", "fracMaxTICBg",0.2,0, 1, 0.05)),
            shiny::column(width = 4,shiny::numericInput("derivThresholdExp", "derivThresholdExp", 1, 0, 5, 0.05)),
            shiny::column(width = 4,shiny::numericInput("derivThresholdBg", "derivThresholdBg",0.05,0, 2, 0.05)),
            shiny::column(width = 4,shiny::numericInput("minPoints", "minPoints", 1, 1, 10, 1)),
            shiny::column(width = 4,shiny::numericInput("degreeBaseline", "degreeBaseline", 1, 1, 30, 1))
        ),
        shiny::fluidRow(
          shiny::p(class='text-center', shiny::actionButton('delete', 'Delete selected rows'),
                   shiny::actionButton('changeTimeLimit', 'Change time limits'))
        )
      )
    }
    
  })
  
  shiny::observeEvent(input$fileNameTimeLimit, {
    rv$dataTimeLimit <- t(rv$ptrset@timeLimit[[input$fileNameTimeLimit]]$exp)
  })
  
  shiny::observeEvent(input$changeTimeLimit, {
    if(!is.null(rv$ptrset) & !is.null(input$fileNameTimeLimit)){
    #shiny::removeNotification("warnings")
      ptrSet<-rv$ptrset
    fracMaxTICNew <- input$fracMaxTIC
    ##default value
    fracMaxTICBgNew <- input$fracMaxTICBg
    derivThresholdExpNew <- input$derivThresholdExp
    derivThresholdBgNew <-input$derivThresholdBg
    minPointsNew <-input$minPoints
    degreeBaselineNew <- input$degreeBaseline
    
    timeLimit<-ptairMS:::timeLimitFun(TIC = ptrSet@breathTracer[[input$fileNameTimeLimit]],
                                      fracMaxTIC =fracMaxTICNew ,fracMaxTICBg = fracMaxTICBgNew,
                                      derivThresholdExp = derivThresholdExpNew,derivThresholdBg = derivThresholdBgNew,
                                      minPoints =minPointsNew ,degreeBaseline = degreeBaselineNew)
    rv$dataTimeLimit <- t(timeLimit$exp)
    rv$bgPoints<- timeLimit$backGround
    #ptrSet<-ptrSetNew()
    ptrSet@timeLimit[[input$fileNameTimeLimit]]$exp <- t(rv$dataTimeLimit)
    ptrSet@timeLimit[[input$fileNameTimeLimit]]$backGround <- rv$bgPoints
    paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                            fracMaxTICBg=fracMaxTICBgNew, 
                            derivThresholdExp=derivThresholdExpNew,
                            derivThresholdBg=derivThresholdBgNew,
                            minPoints = minPointsNew,
                            degreeBaseline=degreeBaselineNew)
    
    ptrSet@parameter$timeLimit<-paramterTimeLimit
    rv$ptrset<-ptrSet
    #ptrSetNew(ptrSet)
    }
  })
  
  output$tableTimeLimit <- DT::renderDataTable( {
    if(!is.null(rv$ptrset) & !is.null(input$fileNameTimeLimit)){
      #fileName<-input$fileNameTimeLimit
      #indexTimeLimit <- rv$ptrset@timeLimit[[fileName]]$exp
      DT::datatable(rv$dataTimeLimit, editable = TRUE,
                    selection=list(target="row"),
                    rownames=as.character(seq(1,nrow(rv$dataTimeLimit)))
      )
    }
  })
  
  
  shiny::observeEvent(input$tableTimeLimit_cell_edit,
                      {
                        rv$dataTimeLimit[input$tableTimeLimit_cell_edit$row, input$tableTimeLimit_cell_edit$col] <- 
                          input$tableTimeLimit_cell_edit$value
                        #ptrSet<-ptrSetNew()
                        rv$ptrset@timeLimit[[input$fileNameTimeLimit]]$exp <- t(rv$dataTimeLimit)
                        #ptrSetNew(ptrSet)
                      }
  )
  
  output$plotTimeLimit <- plotly::renderPlotly({
    if(!is.null(rv$ptrset) & !is.null(input$fileNameTimeLimit)){
      #s <- input$table_rows_selected
      fileName <- input$fileNameTimeLimit
      TIC <- rv$ptrset@breathTracer[[fileName]]
      time <- as.numeric(names(TIC))
      #indexTimeLimit <- t(rv$data)
      indexTimeLimit <- rv$ptrset@timeLimit[[fileName]]$exp
      expPoint<-Reduce(c,apply(indexTimeLimit,2,function(x) seq(x[1],x[2])))
      #bgPoint<-rv$bgPoints
      bgPoint<-rv$ptrset@timeLimit[[fileName]]$backGround
      p <- ggplot2::qplot(x=time,y=TIC,
                          xlab="time",ylab="intensity",main=paste("Mz BreathTracer of",fileName))
      if(!is.null(expPoint))
        p<-p+ggplot2::geom_point(mapping=ggplot2::aes(time ,y,color=point),
                                 data=data.frame(time= time[expPoint],y=TIC[expPoint],point="exp"))
      if(!is.null(bgPoint)) 
        p<-p + ggplot2::geom_point(mapping=ggplot2::aes(time ,y,color=point),
                                   data=data.frame(time= time[bgPoint],y=TIC[bgPoint],point="Background"))
      # if( length(s) ){
      #   s <- s[s<=ncol(indexTimeLimit)] #to avoid warnings
      #   if(length(s)) p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = time[c(indexTimeLimit[,s])],colour="selected"))
      # }
      plotly::ggplotly(p)
    }
   
  })
  
  output$paramCalib<-renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        shiny::fluidRow(
          shiny::selectInput("fileNameCalib","File name:",as.list(names(rv$ptrset@TIC)),selectize=FALSE)
        ),
        shiny::fluidRow(shiny::h4("Calibration masses"),
                        shiny::column(4,shiny::numericInput(inputId = "mz1New",label = "mz1",value = rv$ptrset@parameter$mzCalibRef[1])),
                        shiny::column(4,shiny::numericInput(inputId = "mz2New",label = "mz2",value =rv$ptrset@parameter$mzCalibRef[2])),
                        shiny::column(4,shiny::numericInput(inputId = "mz3New",label = "mz3",value = rv$ptrset@parameter$mzCalibRef[3])),
                        shiny::column(4,shiny::numericInput(inputId = "mz4New",label = "mz4",value = rv$ptrset@parameter$mzCalibRef[4]))
                        ),
        shiny::fluidRow(
          shiny::p(class='text-center', 
                   shiny::actionButton("resetCalib",label = "reset calibration"))
        )
      )
    }
  })
  
  shiny::observeEvent(input$resetCalib,{
    rv$ptrset<-calibration(rv$ptrset,mzCalibRef = c(input$mz1New,input$mz2New,input$mz3New,input$mz4New))
  })
  
  
  output$plotCalib<- renderPlot({
    if(!is.null(rv$ptrset)) plotCalib(rv$ptrset,file=input$fileNameCalib)
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
     if(input$filterExp){
       bgThreshold<-2
       pValThres<-2e-26
     } else {
       bgThreshold<-0
       pValThres<-1
     }
     
     eSet<- ptairMS::alignSamples(ptrSet,ppmGroup = 90,fracGroup = 1- input$fracGroup ,group = NULL,
                                  bgThreshold = bgThreshold,pValGreaterThres = pValThres,pValLessThres = pValThres,
                                  fracExp = input$fracExp,quantiUnit=input$quantiUnit)
     rv$eset<-eSet
  
   })
   
   output$nbPeak <- renderText({
     if(!is.null( rv$eset)){
     X<-Biobase::exprs(rv$eset)
     paste(nrow(X),"peaks aligned")
     }
   })
   
   output$tableAlign<- plotly::renderPlotly({
     if(!is.null( rv$eset)){
     X<-Biobase::exprs(rv$eset)
     
     if(input$Transformation == "log2"){
       X<-log2(X)
     } else if (input$Transformation == "Centred reduced"){
       X<-t(apply(X,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
     }
     #ropls::imageF(log2(Biobase::exprs(rv$eset)))
     plotly::plot_ly(z = X,
                     x = colnames(Biobase::exprs(rv$eset)),
                     y =  paste0("mz-",as.character(rownames(Biobase::exprs(rv$eset)))), type = "heatmap", colorscale='Jet')
     }
   })
   
   
   shiny::observeEvent(input$impute,{
     if(!is.null(rv$eset) & !is.null(rv$ptrset)){
       rv$eset <- impute(eSet = rv$eset,ptrSet = rv$ptrset,parallelize = T)
         }
   })
                       
   
   output$secondSelectionModeledData<-shiny::renderUI({
     if(!is.null(rv$ptrset)){
       shiny::fluidPage( shiny::selectInput("Filemodeled","File to plot",
                                            c("All",rownames(rv$ptrset@sampleMetadata)),
                                            selectize=FALSE,multiple=TRUE),
                         shiny::selectInput("column","Colored by:",
                                            c("files",colnames(rv$ptrset@sampleMetadata)),
                                            selectize=FALSE) )
     }

   })
   
   output$secondSelectionLongitudinal<-shiny::renderUI({
     if(!is.null(rv$ptrset) & !is.null(rv$eset)){
       shiny::fluidPage( 
         shiny::selectInput("mzLongotudinal","Feature", 
                            rownames(Biobase::exprs(rv$eset)),
                            selectize=FALSE),
         shiny::selectInput("idPatient","idPatient", 
                            c("None",colnames(rv$ptrset@sampleMetadata)),
                            selectize=FALSE),
         shiny::selectInput("timeColumn","Time scale ", 
                            c("None",colnames(rv$ptrset@sampleMetadata)),
                            selectize=FALSE),
         shiny::selectInput("columnLongitudinal","Colored by:",
                                            c("None",colnames(rv$ptrset@sampleMetadata)),
                                            selectize=FALSE) )
     }
     
   })
   
  
   output$dataLongitudinal<-plotly::renderPlotly({
     if(!is.null(rv$eset) ){
       if(!is.null(input$timeColumn) & !is.null(input$columnLongitudinal)){
         if(input$timeColumn != "None" & input$columnLongitudinal != "None"){
           X<-Biobase::exprs(rv$eset)
           Y<-Biobase::pData(rv$eset)
           mz<-input$mzLongotudinal
           colTime<-input$timeColumn
           colorby<-input$columnLongitudinal
           name<-input$idPatient
           Y$file<-rownames(Y)
           quanti<-Biobase::annotation(rv$eset)
           
           g <- ggplot2::ggplot()
           for(j in unique(Y[,name])){
             if(all(!is.na(Y[Y[,name] == j,colTime]))){
               g <- g+ suppressWarnings(ggplot2::geom_line(ggplot2::aes(x=timeC,y=cps,color=filter,label=files,group=name),
                                         data=data.frame(timeC =Y[Y[,name] == j,colTime],
                                                         cps=X[as.character(mz),Y[,name]==j],
                                                         filter=Y[Y[,name] == j,colorby],
                                                         files=rownames(Y)[Y[,name] == j],
                                                         name=j))) +
                 suppressWarnings(ggplot2::geom_point(ggplot2::aes(x=timeC,y=cps,color=filter,label=files),
                                     data=data.frame(timeC=Y[Y[,name] == j,colTime],
                                                     cps=X[as.character(mz),Y[,name]==j],
                                                     filter=Y[Y[,name] == j,colorby],
                                                     files=rownames(Y)[Y[,name] == j],
                                                     name=j)))
             }
            
           }
           
           g<-g + ggplot2::ggtitle(mz,annotateVOC(as.numeric(mz),ppm=100)[3]) +
             ggplot2::labs(y= quanti, x = colTime )
           
           
           plotly::ggplotly(g) %>%
             plotly::layout(title = list(text = paste0(mz,
                                                       '<br>',
                                                       '<sup>',
                                                       annotateVOC(as.numeric(mz),ppm=100)[3],
                                                       '</sup>')))
         }
       }
     }
   })
   
                       
  output$pca<-shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        shiny::selectInput("colorPca","Colored by:",
                           colnames(rv$ptrset@sampleMetadata),
                           selectize=FALSE) 
      )
    }
  })
  
   #plot acps
   output$acp <- plotly::renderPlotly({
     if(!is.null(rv$eset)){
       
       # getting Shiny parameters
       
       eset <- rv$eset

       # log2 transformation
       
       Biobase::exprs(eset) <- log2(Biobase::exprs(eset))
       set.pca <- ropls::opls(eset, predI = 4, info.txtC = "none", fig.pdfC = "none")
       rv$pca<-set.pca
       if( !is.na(input$x) & !is.na(input$y)){
         if(input$x !=input$y){
           score_plotly(set.pca,
                        components.vi = c(input$x, input$y),
                        #label.c = label,
                        color.c = input$colorPca,
                        ellipse.l = FALSE,
                        palette.c = "",
                        #info.vc = c("sampleNames",color),
                        size.ls = list(legend_title.i = 8,
                                       legend_text.i = 8))
         }
       }
       
       
     }
   })
   
   output$getLoading<-renderTable({
     if(!is.null(rv$pca)){
       X<-Biobase::exprs(rv$eset)
       importance<-ropls::getLoadingMN(rv$pca)
       tx<-rownames(X)[order(importance[,input$x]^2,decreasing = TRUE)[1:input$maxRow]]
       ty<-rownames(X)[order(importance[,input$y]^2,decreasing = TRUE)[1:input$maxRow]]
       data<-cbind(mz.tx=tx,annotation.tx=ptairMS::annotateVOC(as.numeric(tx))[,3],
                   mz.ty=ty,annotation.ty=ptairMS::annotateVOC(as.numeric(ty))[,3])
       
       colnames(data)<-gsub(pattern = "x",replacement = input$x,x = colnames(data) )
       colnames(data)<-gsub(pattern = "y",replacement = input$y,x = colnames(data) )
       data
     }
     
   })

  
}

#ptrSet <- shiny::runApp(list(ui=ui,server=server))


# Run the application 
shinyApp(ui = ui, server = server)


