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
                         legend_text.i = 13)
    
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
                                         paste0(", color = ", color.c,",fill=",color.c),
                                         ""),
                                  ifelse(label.c != "",
                                         paste0(", label = ", label.c),
                                         ""),
                                  ", text = text))")))
    
    ## text/points [geom_text/geom_point]
    
    if (label.c != "") {
        p <- p + ggplot2::geom_text(size = size.ls[["label.i"]], ggplot2::aes(fontface = "bold"),show.legend=FALSE)
    } else {
        p <- p + ggplot2::geom_point(size = size.ls[["point.i"]], colour="black",
                                     shape = 21)
    }
    
    ## horizontal and vertical lines [geom_hline, geom_vline]
    
    #p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    # ggplot2::geom_vline(ggplot2::aes(xintercept = 0))
    
    ## ellipses [stat_ellipse]
    
    
    if (ellipse.l && color.c != "" && color_type.c == "qualitative")
        p <- p + eval(parse(text = paste0("ggplot2::stat_ellipse(ggplot2::aes(x = ",
                                          paste0('p', components.vi[1]),
                                          ", y = ",
                                          paste0('p', components.vi[2]),
                                          ", group = ",
                                          ifelse(color.c != "", color.c, 1),
                                          "), type = 'norm',geom = 'polygon',alpha = 0.25)")))
    
    # title and axis labels [labs]
    
    p <- p + ggplot2::labs(title = title.c,
                           x = paste0("PC ", components.vi[1],
                                      " (",
                                      round(ropls.model@modelDF[components.vi[1], "R2X"] * 100),
                                      "%)"),
                           y = paste0("PC ", components.vi[2],
                                      " (",
                                      round(ropls.model@modelDF[components.vi[2], "R2X"] * 100),
                                      "%)"))
    
    # theme [them_bw, theme]
    
    p <- p + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = size.ls[["title.i"]],face = "bold"),
                       axis.title.x = ggplot2::element_text(size = size.ls[["axis_lab.i"]]),
                       axis.title.y = ggplot2::element_text(size = size.ls[["axis_lab.i"]]),
                       axis.text = ggplot2::element_text(size = size.ls[["axis_text.i"]]),
                       legend.text = ggplot2::element_text( size = size.ls[["legend_text.i"]]),
                       legend.position =c(0.99, 0.99), legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       legend.margin = ggplot2::margin(6, 6, 6, 6), 
                       legend.title = ggplot2::element_blank())
    
    # palette [scale_colour_brewer, scale_colour_gradientn]
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