
#' @param object Numeric matrix 
#' @param exp.design xxx
#' @param pal.name xxxx
#' @param subset A `integrer()` vector of index indicating the indices
#' of rows in the dataset to highlight
#' 
#' @return A boxplot
#' 
#' @author Samuel Wieczorek, Anais Courtier, Enora Fremy
#' 
#' @import highcharter
#' @importFrom grDevices colorRampPalette boxplot.stats
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.exclude
#' 
#' @export
#' 
#' @rdname plot_ds
#' 
boxPlotD <- function(object,
                     exp.design,
                     subset = NULL,
                     pal.name = NULL
                     ){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  
  qData <- assay(object)
  legend <- exp.design$Sample.name
  
  print(subset)
  if(missing(exp.design))
    stop("'exp.design' is missing.")
  
  # In order to view subset, one need the 'idcol' info
  if (!is.null(subset))
    stopifnot(!is.null(idcol(object)))
  
  myColors <- SampleColors(exp.design$Condition, pal.name)
  
  
  # Internal function
  add_var_to_series_list <- function(x, 
                                     series_list, 
                                     key_vector, 
                                     value_vector){
    base::stopifnot(length(key_vector) == length(value_vector))
    base::stopifnot(length(series_list) == length(key_vector))
    series_list[[x]][length(series_list[[x]])+1] <- value_vector[x]
    names(series_list[[x]])[length(series_list[[x]])] <- key_vector[x]
    return(series_list[[x]])
  }
  
  
  # From highcharter github pages:
  hc_add_series_bwpout <- function(hc, value, by, ...) {
    z = lapply(levels(by), function(x) {
      bpstats = boxplot.stats(value[by == x])$stats
      outliers = c()
      for (y in na.exclude(value[by == x])) {
        if ((y < bpstats[1]) | (y > bpstats[5]))
          outliers = c(outliers, list(which(levels(by)==x)-1, y))
      }
      outliers
    })
    hc %>%
      hc_add_series(data = z, type="scatter", ...)
  }
  
  
  gen_key_vector <- function(variable, num_times)
    return(rep(variable, num_times))
   
  
  
  gen_boxplot_series_from_df <- function(value, by, ...){
    value <- base::as.numeric(value)
    by <- base::as.factor(by)
    box_names<- levels(by)
    z = lapply(box_names, function(x) {
      boxplot.stats(value[by==x])$stats
    })
    tmp<- lapply(seq_along(z), function(x){
      var_name_list<- list(box_names[x])
      #tmp0<- list(names(df)[x])
      names(var_name_list)<- "name"
      index<- x-1
      tmp<- list(c(index,  z[[x]]))
      tmp<- list(tmp)
      names(tmp)<- "data"
      tmp_out<- c(var_name_list, tmp)
      
      return(tmp_out)
    })
    
    return(tmp)
  }
  
  
  # Usage: 
  #series<- gen_boxplot_series_from_df(value = df$total_value, by=df$asset_class)
  
  
  ## Boxplot function:
  make_boxplot<- function(value, 
                          by, 
                          chart_title = "Boxplots",
                          chart_x_axis_label = "Values", 
                          show_outliers = FALSE,
                          boxcolors = NULL, 
                          box_line_colors = NULL){
    by <- as.factor(by)
    box_names_to_use <- levels(by)
    series <- gen_boxplot_series_from_df(value = value, by=by)
    
    cols <- boxcolors
    cols2 <- box_line_colors

    # Injecting value 'fillColor' into series list
    kv<- gen_key_vector(variable = "fillColor", length(series)) 
    series2<- lapply(seq_along(series), function(x){ add_var_to_series_list(x = x, series_list = series, key_vector = kv, value_vector = cols) })
    
    if(show_outliers == TRUE){
      hc<- highcharter::highchart() %>%
        highcharter::hc_chart(type = "boxplot", 
                              inverted = FALSE) %>%
        highcharter::hc_title(text = chart_title) %>%
        highcharter::hc_legend(enabled = FALSE) %>%
        highcharter::hc_xAxis(type = "category", 
                              categories = box_names_to_use, 
                              title = list(text=chart_x_axis_label)) %>%
        highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
        highcharter::hc_add_series_list(series2) %>%
        hc_add_series_bwpout(value = value, 
                             by = by, 
                             name = "Outliers", 
                             colorByPoint = TRUE) %>%
        hc_plotOptions(series = list(
          marker = list(
            symbol = "circle"
          ),
          grouping = FALSE
        )) %>%
        highcharter::hc_colors(cols2) %>%
        highcharter::hc_exporting(enabled = TRUE)
      
    } else{
      hc<- highcharter::highchart() %>%
        highcharter::hc_chart(type="boxplot", inverted=FALSE) %>%
        highcharter::hc_title(text=chart_title) %>%
        highcharter::hc_legend(enabled=FALSE) %>%
        highcharter::hc_xAxis(type="category", 
                              categories = box_names_to_use, 
                              label = list(text = chart_x_axis_label)) %>%
        highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
        highcharter::hc_add_series_list(series2) %>%
        hc_plotOptions(series = list(
          marker = list(
            symbol = "circle"
          ),
          grouping=FALSE
        )) %>%
        highcharter::hc_colors(cols2) %>%
        highcharter::hc_exporting(enabled = TRUE)
    }
    hc
  }
  
  
  df <- data.frame(cbind(categ =rep(colnames(qData),nrow(qData)),
                         value = as.vector(apply(qData, 1, function(x) as.vector(x)))))
  df$value<- base::as.numeric(df$value)
  hc <- make_boxplot(value = df$value, 
                    by = df$categ, 
                    chart_title = "",
                    chart_x_axis_label = "Samples",
                    show_outliers = TRUE, 
                    boxcolors = myColors,
                    box_line_colors = "black")
  
  
  
  # Display of rows to highlight (index of row in subset) 
 if(!is.null(subset)){
   idVector <- rowData(object)[ ,idcol(object)]
   pal.tracker <- ExtendPalette(length(subset), "Dark2")
   n = 0
   for(i in subset){
     n = n+1
     dfSubset <- data.frame(y = as.vector(qData[i,],mode='numeric'), 
                            x = as.numeric(factor(names(qData[i,])))-1, 
                            stringsAsFactors = FALSE)
     hc <- hc %>%
       highcharter::hc_add_series(type = "line",
                                  data = dfSubset,
                                  color = pal.tracker[n], 
                                  dashStyle = "shortdot",
                                  name = idVector[i],
                                  tooltip = list(enabled =TRUE,
                                               headerFormat = "",
                                               pointFormat = "{point.series.name} : {point.y: .2f} ") )
     
   }
 }   
  
  hc
  
}