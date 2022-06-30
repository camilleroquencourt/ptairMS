utils::globalVariables(c("time", "intensity", "sd", "filter","cps","x","z",
                         "intensities","timeC","files"))

plotPeakRaw<-function(set,file, mzRange, peaks, temporalEstim,fctFit,
                      ppm = 2000, 
                      palette = c("heat")[1]) {
  
  
  if (length(mzRange) == 1) 
    mzRange <- mzRange + c(-1, 1) * mzRange * ppm * 1e-06
  mz <- rhdf5::h5read(file, name = "FullSpectra/MassAxis")
  mzRange[1] <- max(min(mz), mzRange[1], na.rm = TRUE)
  mzRange[2] <- min(max(mz), mzRange[2], na.rm = TRUE)
  # get the index of the aroud mzRange
  index <- which(mzRange[1] < mz & mz < mzRange[2])
  # object data
  object <- rhdf5::h5read(file, name = "/FullSpectra/TofData", 
                          index = list(index, NULL, NULL, NULL))
  
  
  # time axis
  timeVn <- as.numeric(names(ptairMS:::getTimeInfo(set)$TIC[[basename(file)]]))
  
  # calibrate mass axis
  FirstcalibCoef <- try(rhdf5::h5read(file, "FullSpectra/MassCalibration", 
                                      index = list(NULL,1)))
  attributCalib <- try(rhdf5::h5readAttributes(file, "/FullSpectra"))
  if (!is.null(attr(FirstcalibCoef, "condition")) & is.null(attr(attributCalib, "condition"))) {
    FirstcalibCoef <- matrix(c(attributCalib$`MassCalibration a`, attributCalib$`MassCalibration b`), 
                             ncol = 1)
  }
  
  tof <- sqrt(mz) * FirstcalibCoef[1, 1] + FirstcalibCoef[2, 1]
  # tof<- seq(0,length(mz))
  coefCalib <- ptairMS:::getCalibrationInfo(set)$coefCalib[[basename(file)]][[1]]
  mzVn <- ((tof - coefCalib["b", ])/coefCalib["a", ])^2
  mzVn <- mzVn[index]
  
  # formate object matix
  rawSubMN <- matrix(object, nrow = dim(object)[1], 
                     ncol = prod(utils::tail(dim(object), 2)))
  rawSubMN <- rawSubMN[,seq_len(length(timeVn)),drop=FALSE]
  
  imageMN<-rawSubMN
  dimnames(imageMN) <- list(mzVn, timeVn)
  
  
  ## chr: Current
  bin<-diff(timeVn)[1]
  chromVn <- apply(rawSubMN, 2, function(intVn) sum(intVn, na.rm = TRUE)/bin )
  
  temporal <- ggplot2::ggplot() + 
    ggplot2::geom_point(data=data.frame(time=timeVn,cps=chromVn),
                        mapping = ggplot2::aes(x=time,y=cps)) + ggplot2::theme_classic() 
  temporal<-temporal + ggplot2::geom_line(data=data.frame(time=timeVn,cps=temporalEstim),
                                          mapping = ggplot2::aes(x=time,y=cps))
  
  
  ## sca: Color scale
  df <- data.frame(x = 1, y = seq_len(236), z = seq_len(236))
  sca <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = z), size = 5) + 
    ggplot2::labs(x = NULL, y = NULL)+
    ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_fill_gradientn(colours = rev(grDevices::rainbow(100, end = 4/6)))
  
  
  ## ima: Image
  
  rawM<-ggplot2::ggplot()
  
  data<-data.frame(time=rep(timeVn,each=nrow(imageMN)),
                   mz=rep(as.numeric(rownames(imageMN)),ncol(imageMN)),
                   intensity=c(imageMN))
  
  rawM <- ggplot2::ggplot() + 
    ggplot2::geom_tile(data=data, mapping=ggplot2::aes(time, mz, fill= intensity),show.legend =FALSE) +
    ggplot2::scale_fill_gradientn(colours = rev(grDevices::rainbow(100, end = 4/6)))+
    ggplot2::theme_classic()+
    ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size =12 ),
                   axis.title.y = ggplot2::element_text(size = 12),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text( size = 10),
                   plot.caption = ggplot2::element_text(hjust = 0.5,size = 13) ) 
  
  
  
  ## spe: Spectrum
  timeiNTER<-diff(timeVn)[1]
  specVn <- apply(rawSubMN, 1, function(intVn) mean(intVn, na.rm = TRUE))
  
  spctr<- ggplot2::ggplot()+ 
    ggplot2::geom_point(data=data.frame(mz=as.numeric(rownames(imageMN)),intensities=specVn),
                        mapping = ggplot2::aes(x=mz,y=intensities)) + ggplot2::theme_classic()+  
    ggplot2::coord_flip()
  
  spctr<- spctr + ggplot2::geom_line(data=data.frame(mz=as.numeric(rownames(imageMN)),
                                                     y=peaks( as.numeric(rownames(imageMN)) )),
                                     mapping=ggplot2::aes(x=mz,y=y))
  
  (ggpubr::ggarrange(
    ggpubr::ggarrange(temporal,rawM,nrow=2,ncol=1,align = "v",heights=c(0.25,0.75)),
    ggpubr::ggarrange(sca,spctr,nrow=2,ncol=1,heights=c(0.25,0.75)),
    widths=c(0.75,0.25),ncol=2,nrow=1))
  
}

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
  
  
  # checking arguments and preparing data ----
  
  ## dataset [ExpressionSet] ----
  
  eset <- ropls::getEset(ropls.model)
  
  if (any(dim(eset) < 1))
    stop("ExpressionSet object could not be extracted from the 'ropls.model'")
  
  pdata.df <- Biobase::pData(eset)
  #pdata.df$date<-as.factor(pdata.df$date)
  data.df <- cbind.data.frame(sampleNames = Biobase::sampleNames(eset),
                              pdata.df,
                              stringsAsFactors = FALSE)
  
  ## plotly info ----
  
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
  
  ## score vectors ----
  
  score.mn <- ropls::getScoreMN(ropls.model)
  data.df <- cbind.data.frame(data.df,
                              text = text.vc,
                              score.mn)
  
  ## labels ----
  
  #stopifnot(length(label.c) == 1)
  
  if (label.c != "")
    stopifnot(label.c %in% colnames(data.df))
  
  ## colors ----
  
  #stopifnot(length(color.c) == 1)
  
  if (color.c != "") {
    stopifnot(color.c %in% colnames(data.df))
    if (is.factor(data.df[, color.c]) || is.character(data.df[, color.c])) {
      color_type.c <- "qualitative"
    } else
      color_type.c <- "quantitative"
  }
  
  ## components ----
  
  stopifnot(length(components.vi) == 2)
  stopifnot(max(components.vi) <= ncol(score.mn))
  
  ## sizes ----
  
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
  
  ## text/points [geom_text/geom_point] ----
  
  if (label.c != "") {
    p <- p + ggplot2::geom_text(size = size.ls[["label.i"]], ggplot2::aes(fontface = "bold"))
  } else {
    p <- p + ggplot2::geom_point(size = size.ls[["point.i"]])
  }
  
  ## horizontal and vertical lines [geom_hline, geom_vline] ----
  
  p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0))
  
  ## ellipses [stat_ellipse] ----
  
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
      p <- p + ggplot2::scale_colour_gradientn(colours = rev(grDevices::rainbow(100, end = 4/6)))
  }
  
  # display/saving [plotly::ggplotly, plotly::layout, htmlwidgets::saveWidget, plotly::as_widget]
  
  if (figure.c == "interactive_plotly" || filename_ext.c == "html") {
    
    p <- plotly::ggplotly(p, tooltip = "text")
    
    p <- plotly::layout(p, hoverlabel = list(font = list(size = 20)))
    
    if (filename_ext.c == "html") {
      
      #htmlwidgets::saveWidget(plotly::as_widget(p), figure.c)
      
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

################# ui -----


#' Graphical interface of ptairMS workflow
#' 
#' The whole workflow of ptairMS can be run interactively through a graphical user interface,
#' which provides visualizations (expiration phases, peaks in the raw data, peak table, individual VOCs),
#' quality controls (calibration, resolution, peak shape and evolution of reagent ions depending on time),
#' and exploratory data analysis.
#' @return Shinny app  
#' @examples
#' \dontrun{RunShinnyApp()}
#' @rdname RunShinnyApp
#' @export

RunShinnyApp<-function(){
    
#' @importFrom shinyscreenshot screenshotButton
ui <- shiny::navbarPage("ptairMS",
                        
                        # Application title
                        #titlePanel("ptairMS"),
                        
                        #tabsetPanel(type = "tabs",
                        shiny::tabPanel("Read and check data", 
                                        
                                        shiny::fluidPage(
                                          shiny::h3("Create/Load a PtrSet a from directory"),
                                          shiny::selectInput("createOrLoad", "Choose action", c("loadPtrSet","createPtrSet"),
                                                             multiple=FALSE),
                                          shiny::uiOutput("param")
                                        ),
                                        
                                        
                                        ##UI Affichage du Ptrset ----
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::h3("Your PtrSet:"),
                                            shiny::verbatimTextOutput("showPtrSet",placeholder = FALSE),
                                            shiny::uiOutput("update")
                                          ),
                                          
                                          
                                          ##UI tracé du PtrSet ----
                                          shiny::mainPanel( 
                                            shiny::plotOutput("plotPtrSet"))),
                                        
                                        ##UI Calibration ----
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::h3("Calibration"),
                                            shiny::uiOutput("paramCalib"),
                                            shiny::uiOutput("paramCalib2"),
                                            shiny::uiOutput("paramCalib3"),
                                            shiny::uiOutput("paramCalib4"), 
                                            shiny::uiOutput("addMass_UI"),
                                            # shiny::uiOutput("addMass")
                                            
                                            
                                            
                                            
                                          ),
                                          
                                          shiny::mainPanel( 
                                            shiny::plotOutput("plotCalib"))),
                                        
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::h3("Expirations phases"),
                                            # shiny::uiOutput("nomFile"),
                                            shiny::uiOutput("paramTimeLimit")
                                          ),
                                          shiny::mainPanel(                          
                                            shiny::splitLayout(
                                              DT::dataTableOutput("tableTimeLimit"),
                                              plotly::plotlyOutput("plotTimeLimit")
                                              
                                            ))
                                          
                                        ),
                                        shiny::fluidPage(
                                          shinyscreenshot::screenshotButton()
                                        )
                                        
                        ),
                        
                        shiny::tabPanel("Detect peak",
                                        
                                        shiny::fluidRow(shiny::h3("Detect peak")),
                                        
                                        shiny::fluidRow(
                                          shiny::column( shiny::numericInput("PPM","min ppm sepration",130,0,200,10),width = 3),
                                          shiny::column(shiny::numericInput("minIntensity","min intensity (cps)",5,0,200,10),width = 3),
                                          shiny::column(shiny::numericInput("resMin ","resolution min",2000,0,20000,1000),width = 3),
                                          shiny::column(shiny::numericInput("resMean ","resolution mean",5000,0,20000,1000),width = 3),
                                          shiny::column(shiny::numericInput("resMax","resolution max",10000,0,20000,1000),width = 3)
                                          
                                          
                                        ),
                                        
                                        shiny::fluidRow(
                                          shiny::column( shiny::actionButton("detectPeak", "Detect peaks"),align="center",width = 12)
                                        ),
                                        
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::h3("Single file visualization"),
                                            shiny::uiOutput("secondSelectionFileRawData"),
                                            shiny::uiOutput("secondSelectionMzRawData")
                                          ),
                                          
                                          #Plot
                                          shiny::mainPanel(shiny::plotOutput("RawData"))
                                          
                                        ),
                                        shiny::fluidPage(
                                          shinyscreenshot::screenshotButton()
                                        )
                        ),
                        
                        
                        shiny::tabPanel("Align samples",
                                        # alignement
                                        
                                        
                                        # alignement
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::fluidPage(
                                              shiny::h3("Alignment between samples"),
                                              shiny::selectInput("quantiUnit","Unit:", c("ppb", "ncps", "cps"),selectize=FALSE),
                                              shiny::numericInput("fracGroup","Maximum percentage of missing values per features",value = 0.3,min = 0,max = 1,step = 0.1),
                                              shiny::uiOutput("group"),
                                              shiny::numericInput("fracExp",
                                                                  "Minimum percentage of samples such that features should be significativly 
                                                  diffrent from ambiant air:",value = 0.1,min = 0,max = 1,step = 0.1),
                                              shiny::numericInput("pValThres",
                                                                  "The minimum p-value of the test between ambiant air and expirations:",value = 0.005,min = 0,max = 1,step = 0.1),
                                              shiny::fluidRow(
                                                shiny::column(width = 4, shiny::actionButton('align', 'Align samples')),
                                                shiny::column(width = 4,shiny::verbatimTextOutput("nbPeak",placeholder = TRUE)),
                                                shiny::column(width = 4,shiny::actionButton('export',"Export"))),
                                              shiny::selectInput("Transformation","Matrix transformation:",c("log2","Centred reduced","None"),
                                                                 selected = FALSE),
                                              shiny::fluidRow(
                                                shiny::column(width = 6, shiny::actionButton("impute", label = "Impute")),
                                                shiny::column(width = 6,shiny::actionButton("outlier","Delete outlier")))
                                              
                                            ), width = 5),
                                          
                                          # Show a plot of the generated distribution
                                          shiny::mainPanel( plotly::plotlyOutput("tableAlign"), width = 7)),
                                        
                                        
                                        shiny::sidebarLayout(
                                          shiny::sidebarPanel(
                                            shiny::h3("Metadata"),
                                            shiny::uiOutput("secondSelectionLongitudinal")
                                          ),
                                          
                                          # Show a plot of the generated distribution
                                          shiny::mainPanel(plotly::plotlyOutput("dataLongitudinal"))
                                          
                                        ),
                                        
                                        
                                        shiny::fluidPage(
                                          shinyscreenshot::screenshotButton()
                                        )
                                        
                        ),
                        shiny::tabPanel("Statistical analysis",
                                        shiny::fluidRow(
                                          shiny::sidebarPanel(
                                            shiny::h3("Principal Component Analysis"),
                                            shiny::numericInput("x","Principal component of x-axis",1,1,4,1),
                                            shiny::numericInput("y","Principal component of y-axis",2,1,4,1),
                                            shiny::uiOutput("pca"),
                                            shiny::uiOutput("pcaMarker")
                                          ),
                                          
                                          # Show a plot of the generated distribution
                                          shiny::mainPanel(  plotly::plotlyOutput("acp"))
                                          
                                        ),
                                        shiny::fluidRow(
                                          shiny::sidebarPanel(
                                            shiny::h3("Loading: features order by contribution"),
                                            shiny::numericInput("maxRow","Number of row to display",10,1,100,1)
                                          ),
                                          shiny::mainPanel(shiny::tableOutput("getLoading"))
                                        ),
                                        shiny::fluidPage(
                                          shinyscreenshot::screenshotButton()
                                        )
                                        
                        )
                        #)
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  paramRv<-shiny::reactiveValues(path=NULL)
  
  rv <- shiny::reactiveValues(ptrset=NULL, data=NULL,bgPoints=NULL,new=FALSE,file=NULL,newFile=TRUE)
  
  ## UI création d'un PtrSet ----
  output$param<- shiny::renderUI({
    if(input$createOrLoad == "createPtrSet"){
      rv$new=TRUE
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
          shiny::column(6,shiny::numericInput(inputId = "mzBreathTracer",label = "mzBreathTracer",value = 59.049))
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
  
  ### Choose a directory ----
  shiny::observeEvent(input$selectDir, {
    pth <- utils::choose.dir() #tk_choose.dir
    paramRv$path<-pth
  })
  
  ### Display the choosen directory ----
  output$showDirectorySelect<-shiny::renderPrint({
    if(!is.null(paramRv$path)) print(paramRv$path)
  })
  
  ### Création d'un directory si selection de createPtrSet ----
  shiny::observeEvent(input$createPtrSet, {
    if(!is.null(input$BoolExp)) fracMaxTIC<-0.7 else fracMaxTIC <- 0
    ptrSet<-createPtrSet(dir = paramRv$path , setName = input$ptrName,
                         mzCalibRef = c(input$mz1,input$mz2,input$mz3,input$mz4,input$mz5),
                         fracMaxTIC = fracMaxTIC,
                         mzBreathTracer = input$mzBreathTracer, saveDir =  paramRv$path)
    rv$ptrset<-ptrSet
    rv$new <-TRUE
  })
  
  shiny::observeEvent(input$selectFile, {
    pth <- file.choose()
    load_obj <- function(f){
      env <- new.env()
      nm <- load(f, env)[1]
      env[[nm]]
    }
    rv$ptrset<-load_obj(pth)
    rv$new <- TRUE
    
  })
  
  ## Serveur Affichage du Ptrset ----
  output$showPtrSet<-shiny::renderPrint(
    if(!is.null(rv$ptrset)) print(rv$ptrset)
  )
  
  ## Serveur Tracé du Ptrset ----
  output$plotPtrSet <- shiny::renderPlot({
    if(!is.null(rv$ptrset)) plot(rv$ptrset)
  })
  
  ## Fonction serveur bouton Update Ptrset
  output$update<- shiny::renderUI({
    if(!is.null(rv$ptrset)) shiny::actionButton(inputId = "updatePtrSet",label = "Update ptrSet")
  })
  
  
  shiny::observeEvent(input$updatePtrSet,{
    if(!is.null(rv$ptrset)){
      rv$ptrset<-updatePtrSet(rv$ptrset)
      rv$new <-TRUE
    }
    
  })
  
  
  ##UI output "paramTimeLimit" lorsque on a un ptrset ----
  output$paramTimeLimit<-shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        
        #choose file in the ptrSet
        shiny::fluidRow(
          if (rv$new){
            shiny::selectInput("fileName","File name:",as.list(names(getTimeInfo(rv$ptrset)$TIC)),selected = FALSE,selectize=FALSE)
            
            
          }
          else {
            shiny::selectInput("fileName","File name:",as.list(names(getTimeInfo(rv$ptrset)$TIC)),selected = input$fileName ,selectize=FALSE)
          }
        ),
        
        
        #UI rentrer les valeurs des paramètres
        
        shiny::fluidRow(shiny::h5(shiny::strong("Parameters :"))),
        
        # if (rv$newFile) {
        # shiny::fluidRow(
        # 
        #   shiny::column(width = 6,shiny::numericInput("fracMaxTIC", "FracMaxTIC", 0.8 ,0, 1, 0.05)),
        #   shiny::column(width = 6,shiny::numericInput("fracMaxTICBg", "FracMaxTICBg",0.2,0, 1, 0.05)),
        #   shiny::column(width = 6,shiny::numericInput("derivThresholdExp", "DerivThresholdExp", 1, 0, 5, 0.05)),
        #   shiny::column(width = 6,shiny::numericInput("derivThresholdBg", "DerivThresholdBg",0.05,0, 2, 0.05)),
        #   shiny::column(width = 6,shiny::numericInput("minPoints", "MinPoints", 1, 1, 10, 1)),
        #   shiny::column(width = 6,shiny::numericInput("degreeBaseline", "DegreeBaseline", 1, 1, 30, 1)),
        #   shiny::column(width = 6,shiny::numericInput("knotsPeriod", "KnotsPeriod", 
        #                                               3, 1, 10, 1)),
        #   shiny::column(width = 10,shiny::selectInput(inputId = "methodKnot", 
        #                                             label = "Method for knot location ", 
        #                                             choices = c("Around expiration",
        #                                                         "Uniform"),
        #                                             multiple = FALSE))
        # )
        # }
        
        
        shiny::fluidRow(
          
          
          shiny::column(width = 6,shiny::numericInput("fracMaxTIC", "FracMaxTIC", rv$ptrset@parameter$timeLimit$fracMaxTIC,0, 1, 0.05)),
          shiny::column(width = 6,shiny::numericInput("fracMaxTICBg", "FracMaxTICBg",rv$ptrset@parameter$timeLimit$fracMaxTICBg,0, 1, 0.05)),
          shiny::column(width = 6,shiny::numericInput("derivThresholdExp", "DerivThresholdExp", rv$ptrset@parameter$timeLimit$derivThresholdExp, 0, 5, 0.05)),
          shiny::column(width = 6,shiny::numericInput("derivThresholdBg", "DerivThresholdBg",rv$ptrset@parameter$timeLimit$derivThresholdBg,0, 2, 0.05)),
          shiny::column(width = 6,shiny::numericInput("minPoints", "MinPoints", rv$ptrset@parameter$timeLimit$minPoints, 1, 10, 1)),
          shiny::column(width = 6,shiny::numericInput("degreeBaseline", "DegreeBaseline", rv$ptrset@parameter$timeLimit$degreeBaseline, 1, 30, 1)),
          shiny::column(width = 6,shiny::numericInput("knotsPeriod", "KnotsPeriod",
                                                      rv$ptrset@parameter$knotsPeriod, 1, 10, 1)),
          shiny::column(width = 10,shiny::selectInput(inputId = "methodKnot",
                                                      label = "Method for knot location ",
                                                      choices = c("Around expiration",
                                                                  "Uniform"),
                                                      multiple = FALSE))
          
          
          
        )
        
        ,
        
        # delete des expirations et d'appliquer les changements de paramètres sur les expirations
        shiny::fluidRow(
          
          shiny::actionButton('reset', 'Apply for this file'),
          shiny::actionButton('resetAll', 'Apply for all')),
        shiny::p(class='text-center', shiny::actionButton('delete', 
                                                          'Delete selected rows')
        )
        
      )
    }
    
  },quoted=FALSE)
  
  
  # output$nomFile<-shiny::renderUI({
  #   if(!is.null(rv$ptrset && rv$new) ){
  #     shiny::fluidPage(
  #    
  #       #choose file in the ptrSet
  #       shiny::fluidRow(
  #         shiny::selectInput("fileName","File name:",as.list(names(rv$ptrset@TIC)),selectize=FALSE)
  #       ))
  #   }
  # }
  # )
  

  
  #rv <- shiny::reactiveValues( data = NULL , bgPoints=NULL, ptrset= ptrSet )
  
  # get the time limit of the selected file
  #rv<-list()
  shiny::observeEvent(input$fileName, {
    if(!is.null(rv$ptrset)){
      rv$data <- t(getTimeInfo(rv$ptrset)$timeLimit[[input$fileName]]$exp)
      rv$bgPoints <- getTimeInfo(rv$ptrset)$timeLimit[[input$fileName]]$backGround
      rv$knots<-  getPeaksInfo(rv$ptrset)$knots[[input$fileName]]
      rv$new <-FALSE
      rv$newFile <-TRUE
    }
  })
  
  #display time limit table
  output$tableTimeLimit <- DT::renderDataTable( {
    if(!is.null(rv$ptrset) & !is.null(input$fileName)){
      #fileName<-input$fileNameTimeLimit
      #indexTimeLimit <- rv$ptrset@timeLimit[[fileName]]$exp
      DT::datatable(rv$data, editable = TRUE,
                    selection=list(target="row"),
                    rownames=as.character(seq(1,nrow(rv$data)))
      )
    }
  })
  ### delete the selected expirations ----
  shiny::observeEvent(input$delete, {
    rowNum <- input$tableTimeLimit_rows_selected
    if(length(rowNum) == nrow(rv$data)){
      shiny::showNotification("At least one period must be selected",
                              duration = NULL, id="warnings" )
    }else {
      #shiny::removeNotification("warnings")
      expirationTodelete<-rv$data[rowNum,,drop=FALSE]
      knots<-rv$knots
      knotTodelete<-Reduce(c,apply(expirationTodelete,1,
                                   function(x) which(x[1] <= knots & 
                                                       knots <= x[2])))
      rv$data <- rv$data[-rowNum,,drop=FALSE]
      rv$knots<-knots[-knotTodelete]
      
      ptrSetTemp<-rv$ptrset
      timeLimitInput<-getTimeInfo(ptrSetTemp)$timeLimit[[input$fileName]]
      timeLimitInput$exp <- t(rv$data)
      ptrSetTemp<-setTimeLimits(ptrSetTemp,newTimeLimites =timeLimitInput,index = input$fileName )
      #ptrSetTemp<-setKnots(ptrSetTemp,rv$knots,input$fileName)
      rv$ptrset<-ptrSetTemp
      #ptrSetNew(ptrSet)
    }
    
  })
  
  ### reset the expirations with timeLimit function for the current file ----
  shiny::observeEvent(input$reset, {
    shiny::removeNotification("warnings")
    fracMaxTICNew <- input$fracMaxTIC
    ##default value
    fracMaxTICBgNew <- input$fracMaxTICBg
    derivThresholdExpNew <- input$derivThresholdExp
    derivThresholdBgNew <-input$derivThresholdBg
    minPointsNew <-input$minPoints
    degreeBaselineNew <- input$degreeBaseline
    
    timeLimit<-timeLimitFun(TIC = getTimeInfo(rv$ptrset)$breathTracer[[input$fileName]],
                            fracMaxTIC =fracMaxTICNew ,
                            fracMaxTICBg = fracMaxTICBgNew,
                            derivThresholdExp = derivThresholdExpNew,
                            derivThresholdBg = derivThresholdBgNew,
                            minPoints =minPointsNew ,
                            degreeBaseline = degreeBaselineNew)
    t<-as.numeric(names( getTimeInfo(rv$ptrset)$TIC[[input$fileName]]))
    background<-timeLimit$backGround
    if(input$methodKnot ==  "Around expiration") 
      method<-"expiration" else method<- "uniform"
    knots<-try(defineKnotsFunc(t,background,input$knotsPeriod,method))
    
    rv$data <- t(timeLimit$exp)
    rv$bgPoints<- timeLimit$backGround
    rv$knots<-knots
    

    timeLimitInput<-getTimeInfo(rv$ptrset)$timeLimit[[input$fileName]]
    
    timeLimitInput$exp <- t(rv$data)
    timeLimitInput$backGround <- rv$bgPoints
    
    rv$ptrset<-setTimeLimits(rv$ptrset,newTimeLimites =timeLimitInput,index = input$fileName )
    rv$ptrset<-setKnots(rv$ptrset,rv$knots,input$fileName)
    
    paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                            fracMaxTICBg=fracMaxTICBgNew, 
                            derivThresholdExp=derivThresholdExpNew,
                            derivThresholdBg=derivThresholdBgNew,
                            minPoints = minPointsNew,
                            degreeBaseline=degreeBaselineNew)
    
    parameter<-getParameters(rv$ptrset)
    parameter$knotsPeriod<-input$knotsPeriod
    parameter$timeLimit<-paramterTimeLimit
    rv$ptrset<- setParameters(rv$ptrset,parameter)

    rv$new <-FALSE
    rv$newFile <-FALSE
  })
  
  
  ### reset the expirations with timeLimit function for all file ----
  
  shiny::observeEvent(input$resetAll, {
    
    #shiny::removeNotification("warnings")
    fracMaxTICNew <- input$fracMaxTIC
    
    ##default value
    fracMaxTICBgNew <- input$fracMaxTICBg
    derivThresholdExpNew <- input$derivThresholdExp
    derivThresholdBgNew <-input$derivThresholdBg
    minPointsNew <-input$minPoints
    degreeBaselineNew <- input$degreeBaseline
    
    
    timeLimit<-lapply(names(getTimeInfo(rv$ptrset)$TIC),function(fileName) {
      timeLimitFun(TIC = getTimeInfo(rv$ptrset)$breathTracer[[fileName]],
                   fracMaxTIC =fracMaxTICNew ,
                   fracMaxTICBg = fracMaxTICBgNew,
                   derivThresholdExp = derivThresholdExpNew,
                   derivThresholdBg = derivThresholdBgNew,
                   minPoints =minPointsNew ,degreeBaseline = degreeBaselineNew)
    })
    
    names(timeLimit) <- names(getTimeInfo(rv$ptrset)$TIC)
    knots<-lapply(names(getTimeInfo(rv$ptrset)$timeLimit),function(fileName) {
      
      background<-timeLimit[[fileName]]$backGround
      t<-as.numeric(names( getTimeInfo(rv$ptrset)$TIC[[fileName]]))
      
      if(input$methodKnot == "Around expiration") {
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
    

    
    rv$ptrset<-setTimeLimits(rv$ptrset,timeLimit)
    rv$ptrset<-setKnots(rv$ptrset,knots)
    
    paramterTimeLimit<-list(fracMaxTIC=fracMaxTICNew,
                            fracMaxTICBg=fracMaxTICBgNew, 
                            derivThresholdExp=derivThresholdExpNew,
                            derivThresholdBg=derivThresholdBgNew,
                            minPoints = minPointsNew,
                            degreeBaseline=degreeBaselineNew)
    
    parameter<-getParameters(rv$ptrset)
    parameter$knotsPeriod<-input$knotsPeriod
    parameter$timeLimit<-paramterTimeLimit
    rv$ptrset<- setParameters(rv$ptrset,parameter)

  })
  
  # display the expirations in data table ----
  output$table <- DT::renderDataTable( {
    if(!is.null(rv$data)) 
    #shiny::removeNotification("warnings")
      DT::datatable(rv$data, editable = TRUE,
                    selection=list(target="row"),
                    rownames=as.character(seq(1,nrow(rv$data)))
      )
  })
  
  shiny::observeEvent(input$tableTimeLimit_cell_edit,
                      {
                        rv$data[ input$tableTimeLimit_cell_edit$row , 
                                 input$tableTimeLimit_cell_edit$col] <- as.numeric(input$tableTimeLimit_cell_edit$value)


                        timeLimitInput<-getTimeInfo(rv$ptrset)$timeLimit[[input$fileName]]
                        timeLimitInput$exp <- t(rv$data)
                        
                        rv$ptrset<-setTimeLimits(rv$ptrset,newTimeLimites =timeLimitInput,index = input$fileName )

                      }
  )
  
  # PlotTimeLimit ----
  output$plotTimeLimit <- plotly::renderPlotly({
    if(!is.null(rv$ptrset) & !is.null(input$fileName) & !is.null(rv$data)){
      s <- input$tableTimeLimit_rows_selected
      fileName <- input$fileName
      TIC <- getTimeInfo(rv$ptrset)$breathTracer[[fileName]]
      time <- as.numeric(names(TIC))
      indexTimeLimit <- t(rv$data)
      #indexTimeLimit <- rv$ptrset@timeLimit[[fileName]]$exp
      expPoint<-Reduce(c,apply(indexTimeLimit,2,function(x) seq(x[1],x[2])))
      bgPoint<-rv$bgPoints
      #bgPoint<-rv$ptrset@timeLimit[[fileName]]$backGround
      knots<-rv$knots
      p <- ggplot2::qplot(x=time,y=TIC,
                          xlab="time",ylab="intensity",
                          main=paste("Trace of",input$fileName))
      if(!is.null(expPoint))
        p<-p+ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                                 data=data.frame(time= time[expPoint],y=TIC[expPoint],point="exp"))
      if(!is.null(bgPoint)) 
        p<-p + ggplot2::geom_point(mapping=ggplot2::aes(time ,y,colour=point),
                                   data=data.frame(time= time[bgPoint],y=TIC[bgPoint],point="Background"))
      
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
    }
    
    
  })
  
  # Calibration ----
  output$paramCalib<- shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        shiny::fluidRow(
          shiny::selectInput("fileNameCalib","File name:",as.list(names(rv$ptrset@TIC)),selectize=FALSE),
          shiny::fluidRow(shiny::h4("Calibration masses"),
                          shiny::h6("To remove a mass, put its value to 0. At least two masses are required")                          ),
          
        ),
        
        # shiny::fluidRow(shiny::h5(strong(length(rv$ptrset@parameter$mzCalibRef)))),
        
        # shiny::fluidRow(shiny::h4("Calibration masses"),
        #                 
        #                 shiny::column(4,shiny::numericInput(inputId = "mz1New",label = "mz1",value = rv$ptrset@parameter$mzCalibRef[1])),
        #                 shiny::column(4,shiny::numericInput(inputId = "mz2New",label = "mz2",value =rv$ptrset@parameter$mzCalibRef[2])),
        #                 shiny::column(4,shiny::numericInput(inputId = "mz3New",label = "mz3",value = rv$ptrset@parameter$mzCalibRef[3])),
        #                 shiny::column(4,shiny::numericInput(inputId = "mz4New",label = "mz4",value = rv$ptrset@parameter$mzCalibRef[4]))
        # ),
        
      )
    }
  })
  
  
  
  output$paramCalib2<- shiny::renderUI({
    if(!is.null(rv$ptrset)){
      
      
      lapply(1:length(rv$ptrset@parameter$mzCalibRef), function(i) {
        
        output[[paste0('b', i)]] <- shiny::renderUI({
          shiny::fluidRow(
            #Affichage des nouvelles masses
            shiny::numericInput(inputId = paste0("mz", i,"New"), label = paste("mz", i), value = rv$ptrset@parameter$mzCalibRef[i]),
            
            
            # #Boutons pour effacer les masses
            # shiny::actionButton(inputId = paste0("close", i),label="X")
          )    
        })
        
      })
    }
  })
  
  
  output$paramCalib3<- shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        
        shiny::fluidRow(
          shiny::p(class='text-center', 
                   shiny::actionButton("resetCalib",label = "Reset calibration"))
        )
        
        # shiny::fluidRow(
        #   shiny::p(class='text-center',
        #            shiny::actionButton("addCalib",label = "Add calibration"))
        # )
        
        
      )
    }
  })
  
  output$paramCalib4 <- shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        
        shiny::fluidRow(
          shiny::p(class='text-center', 
                   shiny::actionButton("addMassButton","Add mass"))
        )
      )
    }
  })
  
  # masse <- eventReactive(input$addMassButton,{
  #   if(!is.null(rv$ptrset)){
  #     #shiny::numericInput(inputId = "IdAddedMass", label = paste("mz", length(rv$ptrset@parameter$mzCalibRef)+1), value = 60)
  #     #print(input$IdAddedMass)
  #     #rv$ptrset@parameter$mzCalibRef[[length(rv$ptrset@parameter$mzCalibRef)+1]] <- input$IdAddedMass
  #     rv$ptrset@parameter$mzCalibRef[[length(rv$ptrset@parameter$mzCalibRef)+1]] <- 0
  #     
  #   }
  # })
  # 
  # output$addMass_UI <- renderUI({
  #   masse()
  # })
  
  # Reset Calib 
  shiny::observeEvent(input$resetCalib,{
    mzCalibRef = c()
    for (i in 1:length(rv$ptrset@parameter$mzCalibRef)){
      if (eval(parse(text=paste0("input$mz",i,"New!=0")))){
        mzCalibRef <-c(mzCalibRef,eval(parse(text=paste0("input$mz",i,"New"))))
        
      }
    }
    
    rv$ptrset<-calibration(rv$ptrset,mzCalibRef)
  })
  
  # Add Mass initialisée à 0
  shiny::observeEvent(input$addMassButton,{
    rv$ptrset@parameter$mzCalibRef[length(rv$ptrset@parameter$mzCalibRef)+1] <- 0
  })
  
  # shiny::observeEvent( function(i){
  #   
  # 
  #   input$closei
  #   {
  #     rv$ptrset@parameter$mzCalibRef <- rv$ptrset@parameter$mzCalibRef[-length(rv$ptrset@parameter$mzCalibRef)]
  # }
  # }
  #   )
  
  # # Close mass
  # shiny::observeEvent( input$close3, {
  #   print(length(rv$ptrset@parameter$mzCalibRef))
  #   rv$ptrset@parameter$mzCalibRef <- rv$ptrset@parameter$mzCalibRef[-3]
  #  
  # }
  # 
  # )
  
  
  
  output$plotCalib<- shiny::renderPlot({
    
    if(!is.null(rv$ptrset)) plotCalib(rv$ptrset,file=input$fileNameCalib)
  })
  
  #  shiny::observeEvent(input$close, {
  #    shiny::stopApp(ptrSet)
  #  })
  #  
  
  
  # Détection des pics ----
  shiny::observeEvent(input$detectPeak, {
    ptrSet<-rv$ptrset
    resolution<-c(as.numeric(input$resMin),
                  as.numeric(input$resMean),
                  as.numeric(input$resMax))
    rv$ptrset<-detectPeak(rv$ptrset,ppm = as.numeric(input$PPM),
                          minIntensity =as.numeric(input$minIntensity) ,
                          resolutionRange =NULL )
    rv$ptrset<-rv$ptrset
  })
  
  ## Choix du fichier ----
  output$secondSelectionFileRawData<-shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage( shiny::selectInput(inputId = "plotRaw",label = "File to plot",
                                           choices= c(names(rv$ptrset@peakList)),
                                           selectize=FALSE)
      )
    }
    
  }) 
  
  # Choix du feature ----
  output$secondSelectionMzRawData<-shiny::renderUI({
    
    if(!is.null(rv$ptrset) & !is.null(input$plotRaw)){
      shiny::fluidPage( shiny::selectInput(inputId = "plotRawmz",label = "features",
                                           choices= c(rownames(Biobase::exprs(rv$ptrset@peakList[[input$plotRaw]]))),
                                           selectize=FALSE)
      )
      
    }
  })
  
  # Plot de détection de peaks ----
  
  output$RawData<-shiny::renderPlot({
    if(!is.null(rv$ptrset)){
      if(!is.null(input$plotRawmz)){
        peakList<-Biobase::fData(ptairMS:::getPeakList(rv$ptrset)[[input$plotRaw]])
        
        peaksParam <- peakList[as.character(input$plotRawmz),
                               c("Mz","parameterPeak.delta1", 
                                 "parameterPeak.delta2", 
                                 "parameterPeak.height")]
        fctFit <- rv$ptrset@fctFit[[input$plotRaw]]
        peakShape<-rv$ptrset@peakShape[[input$plotRaw]]
        peaks <- function(mz) eval(parse(text = paste0("ptairMS:::",fctFit)))(peaksParam[1,"Mz"], 
                                                                              peaksParam[1,"parameterPeak.delta1"], 
                                                                              peaksParam[1,"parameterPeak.delta2" ], 
                                                                              peaksParam[1,"parameterPeak.height"], 
                                                                              mz, 
                                                                              peakShape)
        
        temporalEstim<- Biobase::exprs(ptairMS:::getPeakList(rv$ptrset)[[input$plotRaw]])[as.character(input$plotRawmz),]
        
        
        fileFullName<- getFileNames(rv$ptrset,fullNames =TRUE )[ getFileNames(rv$ptrset)==input$plotRaw]
        g<-plotPeakRaw(set = rv$ptrset,mzRange = as.numeric(input$plotRawmz),
                       file=fileFullName,peaks = peaks,fctFit=fctFit,temporalEstim = temporalEstim,
                       ppm = 700)
        print(g)
      }
      
    }
    
    
  })
  
  output$group <-shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(shiny::selectInput(inputId = "group",label = "Select group",
                                          choices= c("NULL",colnames(rv$ptrset@sampleMetadata)),
                                          selectize=FALSE))
    }
  })
  
  #  
  shiny::observeEvent(input$align, {
    ptrSet<-rv$ptrset
    pValThres<-input$pValThres
    group<-NULL
    if(input$group!="NULL") group<-input$group
    eSet<- ptairMS::alignSamples(ptrSet,ppmGroup = 90,fracGroup = 1- input$fracGroup ,group = group,
                                 pValGreaterThres = pValThres,pValLessThres = pValThres,
                                 fracExp = input$fracExp,quantiUnit=input$quantiUnit)
    rv$eset<-eSet
    
  })
  
  
  output$nbPeak <- shiny::renderText({
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
        X<-t(apply(X,1,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))
      }
      Y<-Biobase::pData(rv$eset)
      #ropls::imageF(log2(Biobase::exprs(rv$eset)))
      plotly::plot_ly(z = X,
                      x = substr(colnames(Biobase::exprs(rv$eset)),1,7),
                      y =  paste0("mz-",as.character(rownames(Biobase::exprs(rv$eset)))), 
                      type = "heatmap", colorscale='Jet')
    }
  })
  
  
  shiny::observeEvent(input$impute,{
    if(!is.null(rv$eset) & !is.null(rv$ptrset)){
      rv$eset <- impute(eSet = rv$eset,ptrSet = rv$ptrset,parallelize = FALSE)
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
        shiny::selectInput("timeColumn","x-axis", 
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
          
          `%>%` <- plotly::`%>%`
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
  
  
  output$pcaMarker<-shiny::renderUI({
    if(!is.null(rv$ptrset)){
      shiny::fluidPage(
        shiny::selectInput("markerPca","Marker:",
                           c("sampleNames",colnames(rv$ptrset@sampleMetadata)),
                           selectize=FALSE) 
      )
    }
  })
  
  
  #plot acps
  output$acp <- plotly::renderPlotly({
    if(!is.null(rv$eset) & !is.null(input$markerPca) ){
      
      # getting Shiny parameters
      
      eset <- rv$eset
      
      # log2 transformation
      X<-Biobase::exprs(eset)
      X[X<0]<-min(X[X>0],na.rm=FALSE)
      Biobase::exprs(eset) <- log2(X)
      set.pca <- ropls::opls(eset, predI = 2, info.txtC = "none", fig.pdfC = "none",
                             crossvalI =4)
      rv$pca<-set.pca
      if( !is.na(input$x) & !is.na(input$y)){
        if(input$x !=input$y){
          score_plotly(set.pca,
                       components.vi = c(input$x, input$y),
                       label.c = input$markerPca,
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
  
  output$getLoading<-shiny::renderTable({
    if(!is.null(rv$pca)){
      X<-Biobase::exprs(rv$eset)
      importance<-ropls::getLoadingMN(rv$pca)
      tx<-rownames(X)[order(importance[,input$x]^2,decreasing = TRUE)[seq_len(input$maxRow)]]
      ty<-rownames(X)[order(importance[,input$y]^2,decreasing = TRUE)[seq_len(input$maxRow)]]
      data<-cbind(mz.tx=tx,annotation.tx=ptairMS::annotateVOC(as.numeric(tx))[,3],
                  mz.ty=ty,annotation.ty=ptairMS::annotateVOC(as.numeric(ty))[,3])
      
      colnames(data)<-gsub(pattern = "x",replacement = input$x,x = colnames(data) )
      colnames(data)<-gsub(pattern = "y",replacement = input$y,x = colnames(data) )
      data
    }
    
  })
  
  
}


    shiny::shinyApp(ui = ui, server = server) 

}



