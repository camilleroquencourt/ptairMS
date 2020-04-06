
#' Shinny app to modify expiration limits.
#'  
#' This function run a shiny app, where you can check the automatic expirations/head spaces detection, and modify it.  
#'  
#' @param ptrSet a ptrSet 
#' @return the ptrSet object modified
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/exhaledAir",  package = "ptairData")
#' ptrSet <- createPtrSet(directory,setName="ptrSet",mzCalibRef=c(21.022,59.049),fracMaxTIC=0.8)
#' \dontrun{ptrSet <- changeTimeLimits(ptrSet)}
#' @export
changeTimeLimits<-function(ptrSet){
  
  ui <- shiny::fluidPage(
    
    #title
    shiny::titlePanel("View and modify time limits"),
    
    shiny::fluidRow(shiny::p("Select a file (you can scroll with the up/down keys on the keyboard)")),
    
    #choose file in the ptrSet
    shiny::fluidRow(
      shiny::selectInput("fileName","File name:",as.list(names(ptrSet@TIC)),selectize=FALSE)
    ),
    
    shiny::fluidRow(shiny::p("Select the periods (rows) that you want delete. It must stay at least one period.
                      Press 'Reset time limits' to reestablish all time period in the current file.")),
    
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
      shiny::p(class='text-center', shiny::actionButton('delete', 'Delete selected rows'),
               shiny::actionButton('reset', 'Reset time limits'))
    ),
    
    
    # change parameter of timeLimit
    shiny::fluidRow(
      shiny::numericInput("fracMaxTIC", "FracMaxTIC", ptrSet@parameter$timeLimit$fracMaxTIC,0, 1, 0.05),
               shiny::numericInput("fracMaxTICBg", "fracMaxTICBg",0.2,0, 1, 0.05),
               shiny::numericInput("derivThresholdExp", "derivThresholdExp", 1, 0, 5, 0.05),
               shiny::numericInput("derivThresholdBg", "derivThresholdBg",0.05,0, 2, 0.05),
               shiny::numericInput("minPoints", "minPoints", 1, 1, 10, 1),
               shiny::numericInput("degreeBaseline", "degreeBaseline", 1, 1, 30, 1)
    ),
    
    shiny::fluidRow(shiny::p("When you have done all file, save the changes and the new ptrSet will be return, 
                      and save in saveDir (set parameter of createPtrSet function) if it is not NULL.")),
    # save the ptrSet
    shiny::fluidRow(
      shiny::p(class = 'text-center', shiny::actionButton('download', 'Save changed and exit app'))
    )
    
  )
  
  server <- function(input, output) {
    
    ptrSetNew <- shiny::reactiveVal(value= ptrSet)
    
    rv <- shiny::reactiveValues( data = NULL , bgPoints=NULL)
    
    # get the time limit of the selected file
    shiny::observeEvent(input$fileName, {
      rv$data <- t(ptrSetNew()@timeLimit[[input$fileName]]$exp)
      rv$bgPoints <- ptrSetNew()@timeLimit[[input$fileName]]$backGround
    })
    
    #delete the selected expirations
    shiny::observeEvent(input$delete, {
        rowNum <- input$table_rows_selected
        if(length(rowNum) == nrow(rv$data)){
          shiny::showNotification(warning("At least one period must be selected"),duration = NULL,id="warnings" )
        }else {
          shiny::removeNotification("warnings")
          rv$data <- rv$data[-rowNum,,drop=FALSE]
          ptrSet<-ptrSetNew()
          ptrSet@timeLimit[[input$fileName]]$exp <- t(rv$data)
          ptrSetNew(ptrSet)
        }
  
    })
    
    # reset the expirations with timeLimit function
    shiny::observeEvent(input$reset, {
      shiny::removeNotification("warnings")
      fracMaxTICNew <- input$fracMaxTIC
      ##default value
      fracMaxTICBgNew <- input$fracMaxTICBg
      derivThresholdExpNew <- input$derivThresholdExp
      derivThresholdBgNew <-input$derivThresholdBg
      minPointsNew <-input$minPoints
      degreeBaselineNew <- input$degreeBaseline
      
      timeLimit<-ptairMS:::timeLimitFun(TIC = ptrSet@breathTracer[[input$fileName]],
                                        fracMaxTIC =fracMaxTICNew ,fracMaxTICBg = fracMaxTICBgNew,
                                        derivThresholdExp = derivThresholdExpNew,derivThresholdBg = derivThresholdBgNew,
                                        minPoints =minPointsNew ,degreeBaseline = degreeBaselineNew)
      rv$data <- t(timeLimit$exp)
      rv$bgPoints<- timeLimit$backGround
      ptrSet<-ptrSetNew()
      ptrSet@timeLimit[[input$fileName]]$exp <- t(rv$data)
      ptrSet@timeLimit[[input$fileName]]$backGround <- rv$bgPoints
      paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                              fracMaxTICBg=fracMaxTICBgNew, 
                              derivThresholdExp=derivThresholdExpNew,
                              derivThresholdBg=derivThresholdBgNew,
                              minPoints = minPointsNew,
                              degreeBaseline=degreeBaselineNew)
      
      ptrSet@parameter$timeLimit<-paramterTimeLimit
      ptrSetNew(ptrSet)
    })
    
    # display the expirations in data table
    output$table <- DT::renderDataTable( {
      DT::datatable(rv$data, #editable = TRUE,
                    selection=list(target="row"),
                    rownames=as.character(seq(1,nrow(rv$data)))
      )
    })
    
    #plot the TIC with expirations limits
    output$TIC <- plotly::renderPlotly({
      s <- input$table_rows_selected
      fileName <- input$fileName
      TIC <- ptrSet@breathTracer[[fileName]]
      time <- as.numeric(names(TIC))
      indexTimeLimit <- t(rv$data)
      expPoint<-Reduce(c,apply(indexTimeLimit,2,function(x) seq(x[1],x[2])))
      bgPoint<-rv$bgPoints
      p <- ggplot2::qplot(x=time,y=TIC,
                          xlab="time",ylab="intensity",main=paste("TIC of",input$fileName))  +
        ggplot2::geom_point(mapping=ggplot2::aes(time ,y,color=point),
                             data=data.frame(time= time[expPoint],y=TIC[expPoint],point="exp")) +
        ggplot2::geom_point(mapping=ggplot2::aes(time ,y,color=point),
                             data=data.frame(time= time[bgPoint],y=TIC[bgPoint],point="Background"))
      if( length(s) ){
        s <- s[s<=ncol(indexTimeLimit)] #to avoid warnings
        if(length(s)) p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = time[c(indexTimeLimit[,s])],colour="selected"))
      }
      plotly::ggplotly(p)
    })
    
    # save the ptrSet with change
    shiny::observeEvent(input$download, {
      ptrSet <- ptrSetNew()
      saveDir<-ptrSet@parameter$saveDir
      if(!is.null(saveDir)){
        name<-ptrSet@parameter$name
        changeName <- parse(text=paste0(name,"<- ptrSet "))
        eval(changeName)
        eval(parse(text =  paste0( "save(" ,name ,",file= paste0( saveDir,'/', '",name,".RData '))")))
        message("ptrSet object save in: ", paste0( saveDir,'/', name,".RData"))
      }
      shiny::stopApp(ptrSet)
    })
    
  }
  
  ptrSetNew <- shiny::runApp(list(ui=ui,server=server))
  return(ptrSetNew)
}
