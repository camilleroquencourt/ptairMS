#' Shinny app to modify expiraion
#'
#' @param ptrSet a ptrSet 
#' @examples
#' @export
changeTimeLimit<-function(ptrSet){
  
  ui <- shiny::fluidPage(
    
    #title
    shiny::titlePanel("Time limits"),
    
    #choose file in the ptrSet
    shiny::fluidRow(
      selectInput("fileName","File name:",as.list(names(ptrSet@TIC)))
    ),
    
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
    
    # save the ptrSet
    shiny::fluidRow(
      p(class = 'text-center', shiny::actionButton('download', 'Saved ptrSet'))
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
      rv$data <- rv$data[-rowNum,,drop=FALSE]
      ptrSet<-ptrSetNew()
      ptrSet@timeLimit[[input$fileName]] <- t(rv$data)
      ptrSetNew(ptrSet)
    })
    
    # reset the expirations with timeLimit function
    shiny::observeEvent(input$reset, {
      rv$data <- t(timeLimitFunc(ptrSet@TIC[[input$fileName]],fracMaxTIC = ptrSet@parameter$fracMaxTIC ))
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
      if(length(s)){
        p<-p + geom_vline(aes(xintercept = time[c(indexTimeLimit[,s])],colour="selected"))
      }
      
      plotly::ggplotly(p)
    })
    
    # save the ptrSet with change
    observeEvent(input$download, {
      ptrSet <- ptrSetNew()
      saveDir<-ptrSet@parameter$dir
      name<-ptrSet@parameter$name
      changeName <- parse(text=paste0(name,"<- ptrSet "))
      eval(changeName)
      eval(parse(text =  paste0( "save(" ,name ,",file= paste0( saveDir,'/', '",name,".RData '))")))
    })
    
  }
  
  shiny::shinyApp(ui = ui, server = server)
}




