utils::globalVariables(c("y","point"))
#' Shinny appplication to modify and view expiration limits
#'  
#' This function run a shiny app, where you can check the automatic expiration
#' detection, knots location, and modify it.  
#'  
#' @param ptrSet a ptrSet object
#' @return the ptrSet object modified
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/exhaledAir",  package = "ptairData")
#' ptrSet <- createPtrSet(directory,setName="ptrSet",mzCalibRef=c(21.022,59.049),
#' fracMaxTIC=0.8)
#' \dontrun{ptrSet <- changeTimeLimits(ptrSet)}
#' @export
changeTimeLimits<-function(ptrSet){
  
  ui <- shiny::fluidPage(
    
    #title
    shiny::titlePanel("View and modify expiration time limits"),

    #choose file in the ptrSet
    shiny::fluidRow(
      shiny::selectInput("fileName","File name:",as.list(names(getTimeInfo(ptrSet)$TIC)),selectize=FALSE)
    ),
    
    shiny::fluidRow(shiny::p("Select the periods (rows) that you want delete. 
    It must stay at least one period.
                      You can change parameters and Press 'Reset time limits' 
                             to re estimate time limits.")),
    
    # TIC plot and timeLImit table
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        DT::dataTableOutput("table")
      ),
      shiny::mainPanel( 
        plotly::plotlyOutput("TIC")
      )
    ),
    
    # delete or reset expirations
    shiny::fluidRow(
      shiny::p(class='text-center', shiny::actionButton('delete', 
                                                        'Delete selected rows'),
               shiny::actionButton('reset', 'apply for this file'),
               shiny::actionButton('resetAll', 'apply for all'))
    ),
    
    
    # change parameter of timeLimit
    shiny::fluidRow(
      shiny::column(width = 3,shiny::numericInput("fracMaxTIC", "FracMaxTIC", 
                                                  getParameters(ptrSet)$timeLimit$fracMaxTIC,0, 1, 0.05)),
      shiny::column(width = 3,shiny::numericInput("fracMaxTICBg", 
                                                  "fracMaxTICBg",0.2,0, 1, 0.05)),
      shiny::column(width = 3,shiny::numericInput("derivThresholdExp", 
                                                  "derivThresholdExp", 1, 0, 5, 
                                                  0.05)),
      shiny::column(width = 3,shiny::numericInput("derivThresholdBg", 
                                                  "derivThresholdBg",0.05,0, 2,
                                                  0.05)),
      shiny::column(width = 3,shiny::numericInput("minPoints", "minPoints", 1, 1, 
                                                  10, 1)),
      shiny::column(width = 3,shiny::numericInput("degreeBaseline", 
                                                  "degreeBaseline", 1, 1, 30, 1)),
      shiny::column(width = 3,shiny::numericInput("knotsPeriod", "knotsPeriod", 
                                                  3, 1, 10, 1)),
      shiny::column(width = 3,shiny::selectInput(inputId = "methodKnot", 
                                                 label = "method for knot location ", 
                                                 choices = c("around expiration",
                                                             "uniform"),
                                                 multiple = FALSE))
    ),
    
    shiny::fluidRow(shiny::p("When you have done all file, save the changes and 
    the new ptrSet will be return,and save in saveDir (parameter of createPtrSet 
                             function) if it is not NULL.")),
    # save the ptrSet
    shiny::fluidRow(
      shiny::p(class = 'text-center', shiny::actionButton('download', 
                                                          'Save changed and 
                                                          exit app'))
    )
    
  )
  
  server <- function(input, output) {
    
    ptrSetNew <- shiny::reactiveVal(value= ptrSet)
    
    rv <- shiny::reactiveValues( data = NULL , bgPoints=NULL)
    
    # get the time limit of the selected file
    shiny::observeEvent(input$fileName, {
      rv$data <- t(getTimeInfo(ptrSetNew())$timeLimit[[input$fileName]]$exp)
      rv$bgPoints <- getTimeInfo(ptrSetNew())$timeLimit[[input$fileName]]$backGround
      rv$knots<-  getPeaksInfo(ptrSetNew())$knots[[input$fileName]]
    })
    
    #delete the selected expirations
    shiny::observeEvent(input$delete, {
        rowNum <- input$table_rows_selected
        if(length(rowNum) == nrow(rv$data)){
          shiny::showNotification("At least one period must be selected",
                                  duration = NULL, id="warnings" )
        }else {
          shiny::removeNotification("warnings")
          expirationTodelete<-rv$data[rowNum,,drop=FALSE]
          knots<-rv$knots
          knotTodelete<-Reduce(c,apply(expirationTodelete,1,
                                       function(x) which(x[1] <= knots & 
                                                           knots <= x[2])))
          rv$data <- rv$data[-rowNum,,drop=FALSE]
          rv$knots<-knots[-knotTodelete]
          ptrSet<-ptrSetNew()
          timeLimitInput<-getTimeInfo(ptrSet)$timeLimit[[input$fileName]]
          timeLimitInput$exp <- t(rv$data)
          ptrSet<-setTimeLimits(ptrSet,newTimeLimites =timeLimitInput,index = input$fileName )
          ptrSet<-setKnots(ptrSet,rv$knots,input$fileName)
         
          ptrSetNew(ptrSet)
        }
  
    })
    
    # reset the expirations with timeLimit function for the current file
    shiny::observeEvent(input$reset, {
      #shiny::removeNotification("warnings")
      fracMaxTICNew <- input$fracMaxTIC
      ##default value
      fracMaxTICBgNew <- input$fracMaxTICBg
      derivThresholdExpNew <- input$derivThresholdExp
      derivThresholdBgNew <-input$derivThresholdBg
      minPointsNew <-input$minPoints
      degreeBaselineNew <- input$degreeBaseline
      
      timeLimit<-timeLimitFun(TIC = getTimeInfo(ptrSet)$breathTracer[[input$fileName]],
                              fracMaxTIC =fracMaxTICNew ,
                              fracMaxTICBg = fracMaxTICBgNew,
                              derivThresholdExp = derivThresholdExpNew,
                              derivThresholdBg = derivThresholdBgNew,
                              minPoints =minPointsNew ,
                              degreeBaseline = degreeBaselineNew)
      t<-as.numeric(names( getTimeInfo(ptrSet)$TIC[[input$fileName]]))
      background<-timeLimit$backGround
      if(input$methodKnot ==  "around expiration") 
        method<-"expiration" else method<- "uniform"
      knots<-try(defineKnotsFunc(t,background,input$knotsPeriod,method))

      rv$data <- t(timeLimit$exp)
      rv$bgPoints<- timeLimit$backGround
      rv$knots<-knots
      
      ptrSet<-ptrSetNew()
      timeLimitInput<-getTimeInfo(ptrSet)$timeLimit[[input$fileName]]
      
      timeLimitInput$exp <- t(rv$data)
      timeLimitInput$backGround <- rv$bgPoints
      
      ptrSet<-setTimeLimits(ptrSet,newTimeLimites =timeLimitInput,index = input$fileName )
      ptrSet<-setKnots(ptrSet,rv$knots,input$fileName)
      
      paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                              fracMaxTICBg=fracMaxTICBgNew, 
                              derivThresholdExp=derivThresholdExpNew,
                              derivThresholdBg=derivThresholdBgNew,
                              minPoints = minPointsNew,
                              degreeBaseline=degreeBaselineNew)
      
      parameter<-getParameters(ptrSet)
      parameter$knotsPeriod<-input$knotsPeriod
      parameter$timeLimit<-paramterTimeLimit
      ptrSet<- setParameters(ptrSet,parameter)
      ptrSetNew(ptrSet)
    })
    
    
    # reset the expirations with timeLimit function for all file
    shiny::observeEvent(input$resetAll, {
      
      #shiny::removeNotification("warnings")
      fracMaxTICNew <- input$fracMaxTIC
      
      ##default value
      fracMaxTICBgNew <- input$fracMaxTICBg
      derivThresholdExpNew <- input$derivThresholdExp
      derivThresholdBgNew <-input$derivThresholdBg
      minPointsNew <-input$minPoints
      degreeBaselineNew <- input$degreeBaseline
      
      
      timeLimit<-lapply(names(getTimeInfo(ptrSet)$TIC),function(fileName) {
       timeLimitFun(TIC = getTimeInfo(ptrSet)$breathTracer[[fileName]],
                    fracMaxTIC =fracMaxTICNew ,
                    fracMaxTICBg = fracMaxTICBgNew,
                    derivThresholdExp = derivThresholdExpNew,
                    derivThresholdBg = derivThresholdBgNew,
                    minPoints =minPointsNew ,degreeBaseline = degreeBaselineNew)
        })
      
      names(timeLimit) <- names(getTimeInfo(ptrSet)$TIC)
      knots<-lapply(names(getTimeInfo(ptrSet)$timeLimit),function(fileName) {
        
        background<-timeLimit[[fileName]]$backGround
        t<-as.numeric(names( getTimeInfo(ptrSet)$TIC[[fileName]]))

        if(input$methodKnot == "around expiration") {
          defineKnotsExpiration(t,background,input$knotsPeriod) 
        } else {
          unique(c(t[1],
                           seq(t[1],utils::tail(t,1),by=input$knotsPeriod),
                           utils::tail(t,1)))
        }
       
      })
      names(knots)<-names(timeLimit)
      
      rv$data <- t(timeLimit[[input$fileName]]$exp)
      rv$bgPoints<- timeLimit[[input$fileName]]$backGround
      rv$knots<-knots[[input$fileName]]
      
      ptrSet<-ptrSetNew()
      
      ptrSet<-setTimeLimits(ptrSet,timeLimit)
      ptrSet<-setKnots(ptrSet,knots)
      
      paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                              fracMaxTICBg=fracMaxTICBgNew, 
                              derivThresholdExp=derivThresholdExpNew,
                              derivThresholdBg=derivThresholdBgNew,
                              minPoints = minPointsNew,
                              degreeBaseline=degreeBaselineNew)
      
      parameter<-getParameters(ptrSet)
      parameter$knotsPeriod<-input$knotsPeriod
      parameter$timeLimit<-paramterTimeLimit
      ptrSet<- setParameters(ptrSet,parameter)
      ptrSetNew(ptrSet)
    })
    
    # display the expirations in data table
    output$table <- DT::renderDataTable( {
      if(!is.null(rv$data)) 
        DT::datatable(rv$data, editable = TRUE,
                    selection=list(target="row"),
                    rownames=as.character(seq(1,nrow(rv$data)))
      )
    })
    
    shiny::observeEvent(input$table_cell_edit,
      {
        rv$data[ input$table_cell_edit$row , 
                 input$table_cell_edit$col] <- as.numeric(input$table_cell_edit$value)
        ptrSet <- ptrSetNew()
        timeLimitInput<-getTimeInfo(ptrSet)$timeLimit[[input$fileName]]
        timeLimitInput$exp <- t(rv$data)

        ptrSet<-setTimeLimits(ptrSet,newTimeLimites =timeLimitInput,index = input$fileName )
        ptrSetNew(ptrSet)
      }
    )
    
    #plot the TIC with expiration limits
    output$TIC <- plotly::renderPlotly({
      s <- input$table_rows_selected
      fileName <- input$fileName
      TIC <- getTimeInfo(ptrSet)$breathTracer[[fileName]]
      time <- as.numeric(names(TIC))
      indexTimeLimit <- t(rv$data)
      expPoint<-Reduce(c,apply(indexTimeLimit,2,function(x) seq(x[1],x[2])))
      bgPoint<-rv$bgPoints
      knots<-rv$knots
      p <- ggplot2::qplot(x=time,y=TIC,
                          xlab="time",ylab="intensity",
                          main=paste("Trace of",input$fileName))
      if(!is.null(expPoint))
        p<-p+ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                             data=data.frame(time= time[expPoint],
                                             y=TIC[expPoint],point="Expiration"))
      if(!is.null(bgPoint)) 
        p<-p + ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                             data=data.frame(time= time[bgPoint],
                                             y=TIC[bgPoint],
                                             point="Background"))
      
      p <- p + ggplot2::theme_classic() 
      
      if(!is.null(knots)){
        p<- p + ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                                    data=data.frame(time= knots,
                                                    y=rep(0,length(knots)),
                                                    point="knots"),shape=4) +
          ggplot2::ggtitle(paste("Trace of",input$fileName,length(knots)))
      }
        
       if( length(s) ){
        s <- s[s<=ncol(indexTimeLimit)] #to avoid warnings
        expirationSelect<-Reduce(c,apply(indexTimeLimit[,s,drop=FALSE],2,function(x) seq(x[1],x[2])))
        if(length(s)) p <- p + 
          ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                              data=data.frame(time= time[expirationSelect],
                                              y=TIC[expirationSelect],
                                              point="Expiration selected"))
        
        p<- p  + ggplot2::scale_colour_manual(values = c("#F8766D","#7CAE00",
                                                         "#C77CFF","#00BFC4")) 
       }
      
      plotly::ggplotly(p)
    })
    
    # save the ptrSet with change
    shiny::observeEvent(input$download, {
      ptrSet <- ptrSetNew()
      saveDir<-getParameters(ptrSet)$saveDir
      if(!is.null(saveDir)){
        name<-getParameters(ptrSet)$name
        changeName <- parse(text=paste0(name,"<- ptrSet "))
        eval(changeName)
        eval(parse(text =  paste0( "save(" ,name ,",file= paste0( saveDir,'/', '",name,".RData '))")))
        message("ptrSet object save in: ",  saveDir,'/', name,".RData")
      }
      shiny::stopApp(ptrSet)
    })
    
  }
  
  ptrSetNew <- shiny::runApp(list(ui=ui,server=server))
  return(ptrSetNew)
}



