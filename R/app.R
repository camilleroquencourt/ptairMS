
#'Shinny app to modify expiration.
#'  
#'This function run a shiny app, where you can check the automatic expirations/head spaces detection, and modify it.  
#'  
#'@param ptrSet a ptrSet 
#'@return the ptrSet object modified
#'@examples
#'library(ptairData)
#'directory <- system.file("extdata/exhaledAir",  package = "ptairData")
#'ptrSet <- createPtrSet(directory,setName="ptrSet",mzCalibRef=c(21.022,59.049),fracMaxTIC=0.8)
#' \dontrun{ptrSet <- changeTimeLimits(ptrSet)}
#' @export
changeTimeLimits<-function(ptrSet){
  
  ui <- shiny::fluidPage(
    
    #title
    shiny::titlePanel("View and modify time limits"),
    
    shiny::fluidRow(p("Select a file (you can scroll with the up/down keys on the keyboard)")),
    
    #choose file in the ptrSet
    shiny::fluidRow(
      selectInput("fileName","File name:",as.list(names(ptrSet@TIC)),selectize=FALSE)
    ),
    
    shiny::fluidRow(p("Select the periods (rows) that you want delete. It must stay at least one period.
                      Press 'Reset time limits' to reestablish all time period in the current file.")),
    
    # TIC plot and timeLImit table
    shiny::sidebarLayout(
      sidebarPanel(
        DT::dataTableOutput("table")
      ),
      mainPanel( 
        plotly::plotlyOutput("TIC")
      )
    ),
    
    # delete or reset expirations
    shiny::fluidRow(
      p(class='text-center', shiny::actionButton('delete', 'Delete selected rows'),
        shiny::actionButton('reset', 'Reset time limits'))
    ),
    
    shiny::fluidRow(p("When you have done all file, save the changes and the new ptrSet will be return, 
                      and save in saveDir (set parameter of createPtrSet function) if it is not NULL.")),
    # save the ptrSet
    shiny::fluidRow(
      p(class = 'text-center', shiny::actionButton('download', 'Save changed and exit app'))
    )
    
  )
  
  server <- function(input, output) {
    
    ptrSetNew <- shiny::reactiveVal(value= ptrSet)
    
    rv <- shiny::reactiveValues( data = NULL )
    
    # get the time limit of the selected file
    shiny::observeEvent(input$fileName, {
      rv$data <- t(ptrSetNew()@timeLimit[[input$fileName]])
    })
    
    #delete the selected expirations
    shiny::observeEvent(input$delete, {
        rowNum <- input$table_rows_selected
        if(length(rowNum) == nrow(rv$data)){
          shiny::showNotification(warning("At least one period must be selected"),duration = NULL,id="warnings" )
        }else {
          removeNotification("warnings")
          rv$data <- rv$data[-rowNum,,drop=FALSE]
          ptrSet<-ptrSetNew()
          ptrSet@timeLimit[[input$fileName]] <- t(rv$data)
          ptrSetNew(ptrSet)
        }
  
    })
    
    # reset the expirations with timeLimit function
    shiny::observeEvent(input$reset, {
      removeNotification("warnings")
      rv$data <- t(timeLimitFun(ptrSet@TIC[[input$fileName]],fracMaxTIC = ptrSet@parameter$fracMaxTIC ))
      ptrSet<-ptrSetNew()
      ptrSet@timeLimit[[input$fileName]] <- t(rv$data)
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
      TIC <- ptrSet@TIC[[fileName]]
      time <- as.numeric(names(TIC))
      indexTimeLimit <- t(rv$data)
      p <- ggplot2::qplot(x=time,y=TIC,
                          xlab="time",ylab="intensity",main=paste("TIC of",input$fileName))  +
        geom_vline(aes(xintercept = time[c(indexTimeLimit)]))
      if( length(s) ){
        s <- s[s<=ncol(indexTimeLimit)] #to avoid warnings
        if(length(s)) p <- p + geom_vline(aes(xintercept = time[c(indexTimeLimit[,s])],colour="selected"))
      }
      plotly::ggplotly(p)
    })
    
    # save the ptrSet with change
    observeEvent(input$download, {
      ptrSet <- ptrSetNew()
      saveDir<-ptrSet@parameter$saveDir
      if(!is.null(saveDir)){
        name<-ptrSet@parameter$name
        changeName <- parse(text=paste0(name,"<- ptrSet "))
        eval(changeName)
        eval(parse(text =  paste0( "save(" ,name ,",file= paste0( saveDir,'/', '",name,".RData '))")))
        message("ptrSet object save in: ", paste0( saveDir,'/', name,".RData"))
      }
      stopApp(ptrSet)
    })
    
  }
  
  ptrSetNew <- runApp(list(ui=ui,server=server))
  return(ptrSetNew)
}
