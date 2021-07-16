

#' @title Computes the connex components of a graph, represented by an adjacency matrix.
#' 
#' @description This function computes the connex components of a graph, represented by an adjacency matrix.
#'  
#' @param obj.se xxxx
#' 
#' @return A list of items containing the connex components derived from the matrices given in arguments
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' cc <- ComputeConnectedComposants(Exp1_R25_pept[[2]])
#' 
#' @export
#' 
ComputeConnectedComposants <- function(obj.se){
  if (length(metadata(obj.se)$ll.AdjMat) == 0)
  {
    warning("Any adjacency matrix is present in the dataset.")
    return(obj.se)
  }
   
  res <- list(all = NULL,
              onlyShared = NULL,
              onlySpec = NULL)
  
  ll.X <- GetAdjMat(obj.se)
  
  if (!is.null(ll.X$all))
    res$all <- get.pep.prot.cc(ll.X$all)

  
  if (!is.null(ll.X$onlyShared))
    res$onlyShared <- get.pep.prot.cc(ll.X$onlyShared)

  if (!is.null(ll.X$onlySpec))
    res$onlySpec <- get.pep.prot.cc(ll.X$onlySpec)

  return (res)
}




#' @title Build the list of connected components of an adjacency matrix.
#' 
#' @description Method to build the list of connected components of an adjacency matrix
#' 
#' @param X An adjacency matrix
#' 
#' @return A list of Connected Components (CC)
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2') 
#' X <- GetAdjMat(Exp1_R25_pept[[2]])$all
#' ll <- get.pep.prot.cc(X)
#' 
#' @export
#' 
#' @importFrom Matrix %&% crossprod
#' 
#' @importFrom graph graphAM connComp 
#' 
get.pep.prot.cc <- function(X){
  if (is.null(X)){
    warning("The adjacency matrix is NULL")
    return()}
  
  p <- dim(X)[2] # Nb proteins
  q <- dim(X)[1] # Nb peptides
  
  multprot.cc <- singprot.cc <- multprot.cc.pep <- singprot.cc.pep <- NULL
  A <- B <- g <- NULL
  ### Adjacency matrix construction
  # boolean matrix product
  print("Start computing boolean matrix product")
  A <- Matrix::crossprod(X, boolArith = TRUE)
  # A <- as.matrix(t(X) %*% X) 
  # A <- as.matrix(t(X) %&% X)
  print("End of computing boolean matrix product")
  
  # remove self-connecting edges
  diag(A) <- rep(0,p)
  A <- matrix(as.numeric(A[,]), ncol=p) # goes back to classical matrix format
  colnames(A) <- rownames(A) <- colnames(X) # reset pep and prot names
  
  # proteins with no shared peptides
  SingleProt.CC.id <- which(rowSums(A)==0) 
  if (length(SingleProt.CC.id) > 0){
    ### Peptides from single prot CCs
    singprot.cc <- as.list(names(SingleProt.CC.id))
    singprot.cc.pep <- list()
    for(i in 1:length(singprot.cc)){
      peplist <- which(X[,singprot.cc[[i]]]!=0)
      singprot.cc.pep[[i]] <- names(peplist)
    }
  }
  
  
  if (length(SingleProt.CC.id) < nrow(A)){
    B <- A[-SingleProt.CC.id,-SingleProt.CC.id] # matrix with no 1-prot CC
    
    ### Protein CCs
    multprot.cc <- NULL
    g <- graph::graphAM(B, edgemode='undirected', values=NA)
    multprot.cc <- graph::connComp(as(g, 'graphNEL'))
    
    ### Peptides from multiple prot CCs
    multprot.cc.pep <- list()
    for(i in 1:length(multprot.cc)){
      protlist <- multprot.cc[[i]]
      subX <- as.matrix(X[,protlist])
      peplist <- which(rowSums(subX)!=0)
      multprot.cc.pep[[i]] <- names(peplist)
    }
  }
  
  
  ### Merge results into a single list
  prot.cc <- c(multprot.cc, singprot.cc)
  pep.cc <- c(multprot.cc.pep, singprot.cc.pep)
  global.cc <- list()
  for(i in 1:length(prot.cc)){
    prot <- prot.cc[[i]]
    pep <- pep.cc[[i]]
    tmp <- list(prot,pep)
    names(tmp) <- c("proteins", "peptides")
    global.cc[[i]] <- tmp
  }
  
  ### Clean memory and return result
  rm(A,B,g,multprot.cc,singprot.cc,multprot.cc.pep,singprot.cc.pep,prot.cc,pep.cc)
  gc()
  return(global.cc)
}







#' @title Build graph definition
#' 
#' @description Builds a graph as defined by its nodes and edges.
#' 
#' @param The.CC A cc (a list)
#' 
#' @param X An adjacency matrix as given by the function BuildListAdjacencyMatrices
#' 
#' @return A list of two itens: nodes and edges  
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2') 
#' X <- GetAdjMat(Exp1_R25_pept[[2]])$all
#' ll <- get.pep.prot.cc(X)
#' g <- buildGraph(ll[[1]], X)
#' 
#' @export
#' 
buildGraph <- function(The.CC, X){
  
  subX <- as.matrix(X[The.CC$peptides, The.CC$proteins])
  nb.prot <- length(The.CC$proteins)
  nb.pep <- length(The.CC$peptides)
  nb.pep.shared <- length(which(rowSums(subX)>1))
  nb.pep.spec <- length(which(rowSums(subX)==1))
  nb.total = nb.prot + nb.pep
  edge.list <- as.data.frame(which(subX==1, arr.ind=TRUE))
  
  def.grp <- c(rep("shared.peptide", nb.pep), rep("protein", nb.prot))
  def.grp[which(rowSums(subX)==1)] <- 'spec.peptide'
  
  
  nodes <- data.frame(id = 1:nb.total,
                      group = def.grp,
                      label = c(rownames(subX), colnames(subX)),
                      title =  paste0("<p>", 1:nb.total,"<br>Tooltip !</p>"),
                      size = c(rep(10, nb.pep),rep(20, nb.prot)),
                      stringsAsFactors = FALSE
  )
  edges <- data.frame(from=c(edge.list$row), 
                      to=c(edge.list$col+ nb.pep),
                      stringsAsFactors = FALSE)
  
  return(list(nodes=nodes, edges=edges))
}




#' @title Display a connected component.
#' 
#' @description Displays a connected component as a graph with the package 'visNetwork'.
#' 
#' @param g A connex component (a list)
#' 
#' @param obj xxx
#' 
#' @param prot.tooltip A vector of strings which are the text for the tooltips of the proteins.
#' 
#' @param pept.tooltip A vector of strings which are the text for the tooltips of the peptides.
#' 
#' @return A plot  
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2') 
#' X <- GetAdjMat(Exp1_R25_pept[[2]])$all
#' ll <- get.pep.prot.cc(X)
#' g <- buildGraph(ll[[1]], X)
#' display.CC.visNet(g)
#' 
#' @importFrom visNetwork visNetwork visNodes visGroups visOptions
#' 
#' @export
#' 
display.CC.visNet <- function(g, 
                              obj=NULL,
                              prot.tooltip=NULL, 
                              pept.tooltip=NULL){
  
  col.prot <- "#ECB57C"
  col.spec <- "#5CA3F7"
  col.shared <- "#0EA513"
  
  
  visNetwork::visNetwork(g$nodes, g$edges, width = "100%", height = "100%") %>%
    visNetwork::visNodes(shape = "dot") %>%                        # square for all nodes
    visNetwork::visGroups(groupname = "spec.peptide", color = col.spec) %>%    # darkblue for group "A"
    visNetwork::visGroups(groupname = "shared.peptide", color = col.shared) %>%    # darkblue for group "A"
    visNetwork::visGroups(groupname = "protein", color = col.prot, shape = "dot") %>%
    visNetwork::visOptions(highlightNearest = FALSE) %>%
    #visLegend()
    #visPhysics(stabilization = FALSE)%>%
    visNetwork::visEdges(color = "#A9A9A9",width = 2)
  
  
}



#' @title Jitter plot of CC
#' 
#' @description Jitter plot of CC
#' 
#' @param list.of.cc List of cc such as returned by the function get.pep.prot.cc
#' 
#' @return A plot  
#' 
#' @author Thomas Burger
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2') 
#' X <- GetAdjMat(Exp1_R25_pept[[2]])$all
#' ll <- get.pep.prot.cc(X)
#' plotJitterCC(ll)
#' 
#' @export
#' 
#' @importFrom graphics plot
#' 
plotJitterCC <- function(list.of.cc){
  if (is.null(list.of.cc)){return()}
  
  length(list.of.cc) # number of CCs
  cc.summary <- sapply(list.of.cc, function(x){c(length(x[[1]]),length(x[[2]]))})
  rownames(cc.summary) <- c("Nb_proteins","Nb_peptides")
  colSums(cc.summary) # total amount of pep and prot in each CC
  colnames(cc.summary) <- 1:length(list.of.cc)
  cc.summary
  rowSums(cc.summary) # c(number of prot, number of pep)
  
  
  cc.summary <- as.data.frame(t(jitter(cc.summary)))
  graphics::plot(jitter(cc.summary[,2]),jitter(cc.summary[,1]), type="p", xlab="#peptides in CC", ylab="#proteins in CC")
  
}





#' @title Display a a jitter plot for CC
#' 
#' @description Display a jitter plot for connex component with the package highcharter
#' 
#' @param df A dataframe containing xxx
#' 
#' @param clickFunction A JavaScript code to show an information when a point is clicked
#' 
#' @return A plot  
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' df <- data.frame(x=1:10, y=1:10)
#' plotJitter_hc(df)
#' 
#' @export
#' 
#' @import highcharter
#' 
plotJitter_hc <- function(df, clickFunction=NULL){
  
  xtitle <- "TO DO"
  
  if (is.null(clickFunction)){
    clickFunction <- 
      JS("function(event) {Shiny.onInputChange('eventPointClicked', [this.index]+'_'+ [this.series.name]);}")
  }
  
  i_tooltip <- which(startsWith(colnames(df),"tooltip"))
  txt_tooltip <- NULL
  for (i in i_tooltip){
    t <- txt_tooltip <- paste(txt_tooltip,"<b>",gsub("tooltip_", "", 
                                                     colnames(df)[i], 
                                                     fixed=TRUE), 
                              " </b>: {point.", colnames(df)[i],"} <br> ", 
                              sep="")
  }
  
  h1 <-  highchart() %>%
    hc_add_series(data = df, type = "scatter") %>%
    dapar_hc_chart(zoomType = "xy",chartType="scatter") %>%
    hc_legend(enabled = FALSE) %>%
    hc_yAxis(title = list(text="Nb of proteins ic CC")) %>%
    hc_xAxis(title = list(text = "Nb of peptides ic CC")) %>%
    hc_tooltip(headerFormat= '',pointFormat = txt_tooltip) %>%
    hc_plotOptions(series = list( animation=list(duration = 100),
                                  cursor = "pointer", 
                                  point = list( events = list( 
                                    click = clickFunction ) ) ) ) %>%
    dapar_hc_ExportMenu(filename = "plotCC")
  
  
  return(h1)
}


