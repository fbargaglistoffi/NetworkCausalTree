#' @title
#' Plot of a Network Causal Tree

#' @description
#' Plots the Network Causal Tree results

#' General parameters
#' @param NCT A NCT object (output of the NetworkCausalTree function).
#' @param cov_names (Ordered) covariates names vector.
#' @param output If "detection" only point estimates are reported, if
#' "estimation" both estimated effects and variances are reported.
#' @param margins 4 - dimensional vector setting the figure margins. The first element denotes
#' the bottom margin, the second the left side margin, the third the top margin
#'  and the fourth the right side marhin
#' @param digits Number of digits to round effect labels
#'

#' Edge - specific parameters
#' @param edge_width Edge width.
#' @param edge_label_cex Edge's label cex.
#' @param edge_arrow_size Edge's arrow size.
#' @param edge_color Edges color.
#' @param edge_label_family Edge label family.
#' @param edge_label_color Edges' label color.

#' Vertices - specific parameters
#' @param vertex_size Vertex size.
#' @param vertex_size2 Vertex size (for rectangular vertices).
#' @param vertex_shape Vertex shape.
#' @param vertex_label_font Vertex label font.
#' @param vertex_label_cex Vertex Label cex.
#' @param effect_color_nodes The estimated effect according to which nodes are colored
#' It can be "1000", "0100","1110","1101".
#' @param vertex_color String Vector of two colors identifying the
#' extremes of the color palette employed to depict positive
#' estimates of the effect specified in "effect_color_nodes"
#' (Negative values are represented in white )
#' @param vertex_frame_color Vertex's frame color.
#' @param vertex_label_color Vertex label color.

#' Title - specific parameters
#' @param title Title.
#' @param main_cex Title cex.
#' @param adj adj.
#' @param main_color Title color.
#' @param main_font Title font.
#'
#' Legend - specific parameters
#' @param legend_color Legend color.
#' @param legend_title_color Legend's title color.
#'
#' @return A plot representing the estimated Network causal tree.
#' Each node reports - from the top to the bottom
#'- i) the estimated effect 10000
#'  ii) the estimated effect 0100
#'  iii) the size of the sub popuedgeson
#'  Colors provide an intuition on the strength of the effect specified in "effect_color_nodes".
#'
#' @importFrom data.tree Node
#' @importFrom igraph graph_from_data_frame V E make_empty_graph layout_as_tree "E<-" "V<-"
#' @importFrom graphics plot lines text legend par
#' @importFrom grDevices colorRampPalette
#' @importFrom utils read.table write.table
#' @importFrom graphics plot.new
#'
#' @export

plot_NCT <- function(NCT,
                     cov_names,
                     title,
                     output = "estimation",
                     effect_color_nodes = "1000",
                     vertex_color = c("lightskyblue","dodgerblue2"),
                     vertex_size = 32,
                     vertex_size2 = 30,
                     vertex_shape = "rectangle",
                     vertex_frame_color = "black",
                     vertex_label_color = "black",
                     vertex_label_cex = 1,
                     vertex_label_font = 1,
                     edge_label_cex = 1,
                     edge_label_family = "sans",
                     edge_color = "black",
                     edge_width = 0.3,
                     edge_arrow_size = 0.5,
                     edge_label_color = "black",
                     legend_color = "black",
                     legend_title_color = "black",
                     main_font = 1,
                     main_cex = 1,
                     adj = 1,
                     main_color = "black",
                     margins = c(0.1, 0.5, 1.5, 0.7),
                     digits = 4
                     ){

  if (!("FILTER" %in% names(NCT)) ||
      nrow(NCT) == 1 ||
      all(is.na(NCT$FILTER)) ||
      length(unique(na.omit(NCT$FILTER))) <= 1) {
    
    # safe minimal plot so testthat sees no error
    par(mar = c(1,1,2,1))
    plot.new()
    title(title)
    text(0.5, 0.5, "Single leaf tree\n(No splits)", cex = 1.3)
    
    return(invisible(NCT))
  }
  
  options(warn = -1)

  # Clean causal rules
  NCT$NOBS <- NCT$NOBS_EST + NCT$NOBS_TR

  if (length(cov_names) > 10){
  for (p in 10 : length(cov_names)) {
    NCT$FILTER <- gsub(pattern = paste("X.", p, sep = ""),
                       replacement = cov_names[p],
                       x = as.character(NCT$FILTER) )
  }
  }
  for (p in 1 : 9) {
    NCT$FILTER <- gsub(pattern = paste("X.", p, sep = ""),
                       replacement = cov_names[p],
                       x = as.character(NCT$FILTER) )
  }

  NCT$FILTER <- gsub(pattern = "data_tree", replacement = "",
                   x = as.character(NCT$FILTER))
  NCT$FILTER <- gsub(pattern = "[$]", replacement = "",
                   x = as.character(NCT$FILTER))
  NCT$FILTER <- gsub(pattern = "_bin", replacement ="",
                   x = as.character(NCT$FILTER))
  NCT$FILTER[which(NCT$FILTER!="NA")] <- paste0("NA & ", NCT$FILTER[which(NCT$FILTER!="NA")])

  if (
    nrow(NCT) == 1 ||                             
    all(is.na(NCT$FILTER)) ||                        
    length(unique(NCT$FILTER[!is.na(NCT$FILTER)])) <= 1
  ) {
   
    par(mar = margins)
    plot.new()
    title(title)
    text(0.5, 0.5,
         "Single Leaf Tree\n(No splits detected)",
         cex = 1.2)
    
    return(invisible(NCT))
  }
  

  # Reshape the NCT object as a tree dataset
  tree_data <- data.tree::as.Node(NCT, mode = "table",
                          pathName = "FILTER",
                          pathDelimiter = " & ",
                          colLevels = NULL,
                          na.rm = TRUE)

  # Clean tree branches
edges_tree <- as.matrix(data.tree::ToDataFrameNetwork(tree_data))

  edges_names <- NULL
  for (i in 1 : nrow(edges_tree)) {
    edges_names_i <- utils::tail(strsplit(edges_tree[,2], "/")[[i]], n = 1)
    edges_names_i <- gsub(pattern = "X.", replacement = "", x = edges_names_i)
    edges_names = c(edges_names, edges_names_i)
  }

  grafo_tree <- igraph::graph_from_edgelist(edges_tree, directed = TRUE)

  if (output == "estimation") {
    for (i in seq_len(nrow(NCT))) {
      node_name <- as.character(NCT$FILTER[i])
      vid <- which(names(V(grafo_tree)) == node_name)
      if (length(vid) == 1) {
        V(grafo_tree)$TAU1000[vid] <- NCT$EFF1000_EST[i]
        V(grafo_tree)$SE1000[vid]  <- NCT$SE1000_EST[i]
        V(grafo_tree)$TAU1101[vid] <- NCT$EFF1101_EST[i]
        V(grafo_tree)$SE1101[vid]  <- NCT$SE1101_EST[i]
        V(grafo_tree)$TAU1110[vid] <- NCT$EFF1110_EST[i]
        V(grafo_tree)$SE1110[vid]  <- NCT$SE1110_EST[i]
        V(grafo_tree)$TAU0100[vid] <- NCT$EFF0100_EST[i]
        V(grafo_tree)$SE0100[vid]  <- NCT$SE0100_EST[i]
      }
    }
  } else {
    for (i in seq_len(nrow(NCT))) {
      node_name <- as.character(NCT$FILTER[i])
      vid <- which(names(V(grafo_tree)) == node_name)
      if (length(vid) == 1) {
        V(grafo_tree)$TAU1000[vid] <- NCT$EFF1000_EST[i]
        V(grafo_tree)$TAU1101[vid] <- NCT$EFF1101_EST[i]
        V(grafo_tree)$TAU1110[vid] <- NCT$EFF1110_EST[i]
        V(grafo_tree)$TAU0100[vid] <- NCT$EFF0100_EST[i]
      }
    }
  }

  if (effect_color_nodes == "1000") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1000
  }

  if (effect_color_nodes == "1101") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1101
  }

  if (effect_color_nodes == "0100") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU0100
  }

  if (effect_color_nodes == "1110") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1110
  }

  cols <- grDevices::colorRampPalette(c(vertex_color[1],
                             vertex_color[2]))(4)

  qeff <- as.numeric(stats::quantile(V(grafo_tree)$stat
                    [which(V(grafo_tree)$stat >= 0)]))[2 : 4]


  V(grafo_tree)$color <- "white"
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= 0
                      & V(grafo_tree)$stat < qeff[1])] <- cols[1]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[1]
                      & V(grafo_tree)$stat < qeff[2])] <- cols[2]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[2]
                      & V(grafo_tree)$stat < qeff[3])] <- cols[3]
  V(grafo_tree)$color[which(V(grafo_tree)$stat >= qeff[3])] <- cols[4]


  # Clean edges names
 edges_names_final <- rep(NULL, length(edges_names))
 edgelabel <- vector(mode = "list",
                      length = length(edges_names))

  for (l in 1 : length(edges_names)) {

    if (length(grep(x = edges_names[l], pattern = ">=")) > 0) {
      edgelabel[[l]][1] <- strsplit(edges_names[l], "(?=[>=])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(edges_names[l], "(?=[>=])", perl = TRUE)[[1]][2:
                          (lengths(gregexpr(",", strsplit(edges_names[l],"(?=[>=])", perl = TRUE)))
                          + 1)],collapse = " ")
    }

    if (length(grep(x = edges_names[l], pattern = "<=")) > 0) {
      edgelabel[[l]][1] <- strsplit(edges_names[l], "(?=[<=])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(edges_names[l], "(?=[<=])",perl = TRUE)[[1]]
                                 [2 : (lengths(gregexpr(",", strsplit(edges_names[l],
                                  "(?=[<=])", perl = TRUE))) + 1)], collapse = " ")
    }

    if (length(grep(x = edges_names[l], pattern = "<=")) == 0
        & length(grep(x = edges_names[l], pattern = "<")) > 0) {
      edgelabel[[l]][1] <- strsplit(edges_names[l], "(?=[<])", perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(edges_names[l], "(?=[<])", perl = TRUE)[[1]]
                                 [2 : (lengths(gregexpr(",", strsplit(edges_names[l],
                                 "(?=[<])", perl = TRUE))) + 1)], collapse = " ")
    }

    if (length(grep(x = edges_names[l], pattern=">=")) == 0 &
        length(grep(x = edges_names[l], pattern=">")) > 0) {

      edgelabel[[l]][1] <- strsplit(edges_names[l], "(?=[>])",
                                    perl = TRUE)[[1]][1]
      edgelabel[[l]][2] <- paste(value = strsplit(edges_names[l], "(?=[>])",
                                perl = TRUE)[[1]][2: (lengths(gregexpr(",",
                                strsplit(edges_names[l], "(?=[>])", perl = TRUE))) + 1)],
                                collapse = " ")
    }

  edges_names_final[l] <- paste(edgelabel[[l]][1],
                                gsub("\\s", "",
                                x = edgelabel[[l]][2]),
                                sep = "\n")

  }

edges_names_final <- gsub(pattern = "<\\(1\\)", replacement = "=0", x = as.character(edges_names_final))
edges_names_final <- gsub(pattern = ">=\\(1\\)", replacement = "=1", x = as.character(edges_names_final))
E(grafo_tree)$label <- edges_names_final

  if (output == "estimation") {
    eff1 <- paste(round(NCT$EFF1000_EST, digits), "(", round(NCT$SE1000_EST, digits), ")", sep = "")
    eff2 <- paste(round(NCT$EFF1101_EST, digits), "(", round(NCT$SE1101_EST, digits), ")", sep = "")
    eff3 <- paste(round(NCT$EFF1110_EST, digits), "(", round(NCT$SE1110_EST, digits), ")", sep = "")
    eff4 <- paste(round(NCT$EFF0100_EST, digits), "(", round(NCT$SE0100_EST, digits), ")", sep = "")
  }
  else {
    eff1 <- paste(round(NCT$EFF1000_EST, digits), sep="")
    eff2 <- paste(round(NCT$EFF1101_EST, digits), sep="")
    eff3 <- paste(round(NCT$EFF1110_EST, digits), sep="")
    eff4 <- paste(round(NCT$EFF0100_EST, digits), sep="")
  }

  # First line: direct effect τ(1,0;0,0)
  # Second line: spillover effect τ(0,1;0,0)
  # Third line: number of observations
  V(grafo_tree)$labels <- paste(eff1, eff4, NCT$NOBS, sep = "\n")

  # Plot the tree
  par(mar = margins)

  NCTPLOT <- plot(grafo_tree,
                  layout = layout_as_tree(grafo_tree),
                  edge.label.color = edge_label_color,
                  edge.width = edge_width,
                  edge.label.cex = edge_label_cex,
                  edge.label.family = edge_label_family,
                  vertex.color = V(grafo_tree)$color,
                  vertex.label.dist = 0,
                  vertex.label.color = vertex_label_color,
                  edge.arrow.size = edge_arrow_size,
                  vertex.label = V(grafo_tree)$labels,
                  vertex.label.font = vertex_label_font,
                  vertex.label.cex = vertex_label_cex,
                  vertex.frame.color = vertex_frame_color,
                  edge.color = edge_color,
                  vertex.frame = 2,
                  edge.label = E(grafo_tree)$label,
                  vertex.shape = vertex_shape,
                  vertex.size = vertex_size,
                  vertex.size2 = vertex_size2)

  legend('topleft', text.font = 4, text.col = legend_color,
         xjust = 0, adj = 0.2, legend = c(expression(paste(tau, "(1,0;0,0)",sep = "")),
         expression(paste(tau,"(0,1;0,0)", sep = "")), "N"),
         title = "Legend", title.col  = legend_title_color)

  title(title, cex.main = main_cex, col.main = main_color,
        adj = adj, font.main = main_font)
}

