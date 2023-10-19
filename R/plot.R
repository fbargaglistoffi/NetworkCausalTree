#' @title
#' Plot of a Network Causal Tree

#' @description
#' Plots the Network Causal Tree results

#' General parameters
#' @param NCT A NCT object (output of the NetworkCausalTree function).
#' @param varnames (Ordered) covariates names vector.
#' @param output If "detection" only point estimates are reported, if
#' "estimation" both estimated effects and variances are reported.

#' Edge - specific parameters
#' @param ewidth Edge width.
#' @param elabelcex Edge's label cex.
#' @param edge.arrow.size Edge's arrow size.
#' @param ecolor Edges color.
#' @param elabelfamily Edge label family.
#' @param elabelcolor Edges' label color.

#' Vertices - specific parameters
#' @param vsize Vertex size.
#' @param vsize2 Vertex size (for rectangular vertices).
#' @param vshape Vertex shape.
#' @param vlabelfont Vertex label font.
#' @param vlabelcex Vertex Label cex.
#' @param coloreff The estimated effect according to which nodes are colored
#' It can be "1000", "0100","1110","1101".
#' @param vcolor String Vector of two colors identifying the
#' extremes of the color palette employed to depict positive 
#' estimates of the effect specified in "coloreff" 
#' (Negative values are represented in white )
#' @param vframecolor Vertex's frame color.
#' @param vlabelcolor Vertex label color.

#' Title - specific parameters
#' @param title Title.
#' @param cex.main Title cex.
#' @param adj adj.
#' @param col.main Title color.

#' Legend - specific parameters
#' @param colleg Legend color.
#' @param coltitleleg Legend's title color.
#'
#' @return A plot representing the estimated Network causal tree. 
#' Each node reports - from the top to the bottom
#'- i) the estimated effect 10000
#'  ii) the estimated effect 0100
#'  iii) the size of the sub population
#'  Colors provide an intuition on the strength of the effect specified in "coloreff". 

#'
#' @import data.tree
#' @import igraph
#'
#' @export

plot_NCT <- function(NCT,
                     varnames,
                     title,
                     output = "estimation",
                     coloreff = "1000",
                     vcolor = c("lightskyblue","dodgerblue2"),
                     vsize = 44,
                     vsize2 = 40,
                     vshape = "rectangle",
                     vframecolor = "black",
                     vlabelcolor = "black",
                     vlabelcex = 2,
                     vlabelfont = 1,
                     elabelcex = 1.5,
                     elabelfamily = "sans",
                     ecolor = "black",
                     ewidth = 0.3,
                     edge.arrow.size = 0.5,
                     elabelcolor = "black",  
                     colleg = "black",
                     coltitleleg = "black",
                     font.main = 1,
                     cex.main = 1,
                     adj = 1,
                     col.main = "black"){

  options(warn = -1)
  
  # Clean causal rules
  NCT$NOBS <- NCT$NOBS_EST + NCT$NOBS_TR

  if (length(varnames) > 10){
  for (p in 10 : length(varnames)) {
    NCT$FILTER <- gsub(pattern = paste("X.", p, sep=""),
                     replacement = varnames[p], 
                     x = as.character(NCT$FILTER) )
  }
  }  
  for (p in 1 : 9) {
    NCT$FILTER <- gsub(pattern = paste("X.", p, sep=""),
                     replacement = varnames[p], 
                     x = as.character(NCT$FILTER) )
  }
  
  NCT$FILTER<-gsub(pattern = "data_tree", replacement = "",
                   x=as.character(NCT$FILTER))
  NCT$FILTER<-gsub(pattern = "[$]", replacement = "",
                   x=as.character(NCT$FILTER))
  NCT$FILTER<-gsub(pattern = "_bin", replacement ="",
                   x=as.character(NCT$FILTER))
  NCT$FILTER[which(NCT$FILTER!="NA")] <- paste0("NA & ", NCT$FILTER[which(NCT$FILTER!="NA")])

  # Reshape the NCT object as a tree dataset 
  tree_data <- as.Node(NCT, mode = "table",
                       pathName = "FILTER", 
                       pathDelimiter = " & ", 
                       colLevels = NULL,
                       na.rm = TRUE) 
  
  # Clean tree branches
  lati_tree <- ToDataFrameNetwork(tree_data)
  lati_tree <- as.matrix(lati_tree)
  nomelati <- NULL
  for (i in 1 : nrow(lati_tree)) {
    nomelati_i <- tail(strsplit(lati_tree[,2], "/")[[i]], n = 1)
    nomelati_i <- gsub(pattern = "X.", replacement = "", x = nomelati_i)
    nomelati = c(nomelati, nomelati_i)
  } 
  
  # Generate a tree graph corresponding to the NCT object
  grafo_tree <- graph_from_edgelist(lati_tree, directed = TRUE)
  
  # Add the effects as nodes' attributes
  V(grafo_tree)$TAU1000 <- NCT$EFF1000_EST
  V(grafo_tree)$SE1000 <- NCT$SE1000_EST
  V(grafo_tree)$TAU1101 <- NCT$EFF1101_EST
  V(grafo_tree)$SE1101 <- NCT$SE1101_EST
  V(grafo_tree)$TAU1110 <- NCT$EFF1110_EST
  V(grafo_tree)$SE1110 <- NCT$SE1110_EST
  V(grafo_tree)$TAU0100 <- NCT$EFF0100_EST
  V(grafo_tree)$SE0100 <- NCT$SE0100_EST
  
  # Set the nodes' colors and attach them as an attribute 
  if (coloreff == "1000") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1000 
  }
  
  if (coloreff == "1101") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1101
  } 
  
  if (coloreff == "0100") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU0100
  }
  
  if (coloreff == "1110") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1110
  }
  
  cols <- colorRampPalette(c(vcolor[1],
                             vcolor[2]))(4)
  
  qeff <- as.numeric(quantile(V(grafo_tree)$stat
                              [which(V(grafo_tree)$stat >= 0)]))[2 : 4]
  
  
  V(grafo_tree)$color[which(V(grafo_tree)$stat <  (0)) ]<- "white"
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= 0 
                            & V(grafo_tree)$stat < qeff[1])] <- cols[1]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[1] 
                            & V(grafo_tree)$stat < qeff[2])] <- cols[2]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[2] 
                            & V(grafo_tree)$stat < qeff[3])] <- cols[3]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[3])] <- cols[4]
  
  
  # Clean edges names
  latinomefine <- rep(NULL, length(nomelati))
  edgelabel <- vector(mode = "list", 
                      length = length(nomelati))
  
  for (l in 1 : length(nomelati)) {
    
    if (length(grep(x = nomelati[l], pattern = ">=")) > 0) {
      edgelabel[[l]][1] <- strsplit(nomelati[l], "(?=[>=])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(nomelati[l], "(?=[>=])", perl = TRUE)[[1]][2:
                          (lengths(gregexpr(",", strsplit(nomelati[l],"(?=[>=])", perl = TRUE))) 
                          + 1)],collapse = " ")
    }
    
    if (length(grep(x = nomelati[l], pattern = "<=")) > 0) {
      edgelabel[[l]][1] <- strsplit(nomelati[l], "(?=[<=])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(nomelati[l], "(?=[<=])",perl = TRUE)[[1]]
                                 [2 : (lengths(gregexpr(",", strsplit(nomelati[l],
                                  "(?=[<=])", perl = TRUE))) + 1)], collapse = " ")
    }
    
    if (length(grep(x = nomelati[l], pattern = "<=")) == 0 
        & length(grep(x = nomelati[l], pattern = "<")) > 0) {
      edgelabel[[l]][1] <- strsplit(nomelati[l], "(?=[<])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(nomelati[l], "(?=[<])", perl = TRUE)[[1]]
                                 [2 : (lengths(gregexpr(",", strsplit(nomelati[l],
                                 "(?=[<])", perl = TRUE))) + 1)], collapse = " ")
    }
    
    if (length(grep(x = nomelati[l], pattern=">=")) == 0 &
        length(grep(x = nomelati[l], pattern=">")) > 0) {
      
      edgelabel[[l]][1] <- strsplit(nomelati[l], "(?=[>])",
                                    perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(nomelati[l], "(?=[>])", 
                                perl = TRUE)[[1]][2: (lengths(gregexpr(",", 
                                strsplit(nomelati[l], "(?=[>])", perl = TRUE))) + 1)],
                                collapse = " ")
    }
    
    latinomefine[l] <- paste(edgelabel[[l]][1], 
                             gsub("\\s", "", x = edgelabel[[l]][2]), sep = "\n")
    
  }
  
  latinomefine <- gsub(pattern = "<\\(1\\)", replacement = "=0", x = as.character(latinomefine))
  latinomefine <- gsub(pattern = ">=\\(1\\)", replacement = "=1", x = as.character(latinomefine))
  E(grafo_tree)$label <- latinomefine
  
  # Attach nodes' labels
  eff1 <- paste(round(NCT$EFF1000_EST, 2), "(", round(NCT$SE1000_EST, 2), ")", sep="")
  eff2 <- paste(round(NCT$EFF1101_EST, 2), "(", round(NCT$SE1101_EST, 2), ")", sep="")
  eff3 <- paste(round(NCT$EFF1110_EST, 2), "(", round(NCT$SE1110_EST, 2), ")", sep="")
  eff4 <- paste(round(NCT$EFF0100_EST, 2), "(", round(NCT$SE0100_EST, 2), ")", sep="") 
  V(grafo_tree)$labels <- paste(eff1, eff4, NCT$NOBS, sep="\n") 
  
  # Plot the tree
  NCTPLOT <- plot(grafo_tree,
                  layout = layout_as_tree(grafo_tree), 
                  edge.label.color = elabelcolor,
                  edge.width = ewidth, 
                  edge.label.cex = elabelcex, 
                  edge.label.family = elabelfamily,
                  vertex.color = V(grafo_tree)$color,
                  vertex.label.dist = 0,
                  vertex.label.color = vlabelcolor,
                  edge.arrow.size = edge.arrow.size,
                  vertex.label = V(grafo_tree)$labels,  
                  vertex.label.font = vlabelfont, 
                  vertex.label.cex = vlabelcex,
                  vertex.frame.color = vframecolor,
                  edge.color = ecolor, 
                  vertex.frame=2,
                  edge.label = E(grafo_tree)$label,
                  vertex.shape = vshape,
                  vertex.size = vsize,
                  vertex.size2 = vsize2)
  
  legend('topleft', text.font = 4, text.col = colleg, 
         xjust = 0, adj = 0.2, legend = c(expression(paste(tau, "(1,0;0,0)",sep = "")),
         expression(paste(tau,"(0,1;0,0)", sep = "")), "N"),
         title = "Legend", title.col  = coltitleleg)  

  title(title, cex.main = cex.main, col.main = col.main,
        adj = adj, font.main = font.main)
}

