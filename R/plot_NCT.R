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
#'  and the fourth the right side margin
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
#'- i) the estimated effect 1000
#'  ii) the estimated effect 0100
#'  iii) the size of the sub population
#'  Colors provide an intuition on the strength of the effect specified in "effect_color_nodes".
#'
#' @importFrom data.tree Node
#' @importFrom igraph graph_from_data_frame V E make_empty_graph layout_as_tree "E<-" "V<-"
#' @importFrom graphics plot lines text legend par
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot.new
#'
#' @export

plot_NCT <- function(NCT,
                     cov_names,
                     title,
                     output = "estimation",
                     effect_color_nodes = "1000",
                     vertex_color = NULL,
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
                     margins = c(0.1, 1.5, 1.5, 1.5),
                     digits = 4
){
  
  if (!("FILTER" %in% names(NCT)) ||
      nrow(NCT) == 1 ||
      all(is.na(NCT$FILTER)) ||
      length(unique(na.omit(NCT$FILTER))) <= 1) {
    
    par(mar = c(1,1,2,1))
    plot.new()
    title(title)
    text(0.5, 0.5, "Single leaf tree\n(No splits)", cex = 1.3)
    
    return(invisible(NCT))
  }
  
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn), add = TRUE)
  
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
  
  tree_data <- data.tree::as.Node(NCT, mode = "table",
                                  pathName = "FILTER",
                                  pathDelimiter = " & ",
                                  colLevels = NULL,
                                  na.rm = TRUE)
  
  edges_tree <- as.matrix(data.tree::ToDataFrameNetwork(tree_data))
  
  if (nrow(edges_tree) == 0) {
    par(mar = c(1,1,2,1))
    plot.new()
    title(title)
    text(0.5, 0.5, "Single leaf tree\n(No edges)", cex = 1.3)
    return(invisible(NCT))
  }
  
  edges_names <- NULL
  for (i in 1 : nrow(edges_tree)) {
    edges_names_i <- utils::tail(strsplit(edges_tree[,2], "/")[[i]], n = 1)
    edges_names_i <- gsub(pattern = "X.", replacement = "", x = edges_names_i)
    edges_names = c(edges_names, edges_names_i)
  }
  
  grafo_tree <- igraph::graph_from_edgelist(edges_tree, directed = TRUE)
  
  if (output == "estimation") {
    V(grafo_tree)$TAU1000 <- NCT$EFF1000_EST
    V(grafo_tree)$SE1000 <- NCT$SE1000_EST
    V(grafo_tree)$TAU1101 <- NCT$EFF1101_EST
    V(grafo_tree)$SE1101 <- NCT$SE1101_EST
    V(grafo_tree)$TAU1110 <- NCT$EFF1110_EST
    V(grafo_tree)$SE1110 <- NCT$SE1110_EST
    V(grafo_tree)$TAU0100 <- NCT$EFF0100_EST
    V(grafo_tree)$SE0100 <- NCT$SE0100_EST    
  } else {
    V(grafo_tree)$TAU1000 <- NCT$EFF1000_EST
    V(grafo_tree)$TAU1101 <- NCT$EFF1101_EST
    V(grafo_tree)$TAU1110 <- NCT$EFF1110_EST
    V(grafo_tree)$TAU0100 <- NCT$EFF0100_EST
  }
  
  if (effect_color_nodes == "1000") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1000
  } else if (effect_color_nodes == "1101") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1101
  } else if (effect_color_nodes == "0100") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU0100
  } else if (effect_color_nodes == "1110") {
    V(grafo_tree)$stat <- V(grafo_tree)$TAU1110
  }
  
  # Color palette
  if (is.null(vertex_color)) {
    if (effect_color_nodes == "1000" || effect_color_nodes == "1101" || effect_color_nodes == "1110") {
      vertex_color <- c("lightblue", "dodgerblue3")
    } else {
      vertex_color <- c("lightgreen", "forestgreen")
    }
  }
  
  stat_values <- V(grafo_tree)$stat
  positive_stats <- stat_values[stat_values >= 0]
  
  V(grafo_tree)$color <- "white"
  
  if (length(positive_stats) > 0) {
    n_colors <- 15
    cols <- grDevices::colorRampPalette(c(vertex_color[1], vertex_color[2]))(n_colors)
    
    max_stat <- max(positive_stats)
    min_stat <- min(positive_stats)
    
    if (max_stat > min_stat) {
      normalized <- (positive_stats - min_stat) / (max_stat - min_stat)
      color_indices <- pmax(1, pmin(n_colors, ceiling(normalized * n_colors)))
      V(grafo_tree)$color[stat_values >= 0] <- cols[color_indices]
    } else {
      V(grafo_tree)$color[stat_values >= 0] <- cols[ceiling(n_colors/2)]
    }
  }
  
  edges_names_final <- rep(NULL, length(edges_names))
  edgelabel <- vector(mode = "list", length = length(edges_names))
  
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
    eff1 <- paste(round(NCT$EFF1000_EST, digits), "(", round(NCT$SE1000_EST, digits), ")", sep="")
    eff2 <- paste(round(NCT$EFF1101_EST, digits), "(", round(NCT$SE1101_EST, digits), ")", sep="")
    eff3 <- paste(round(NCT$EFF1110_EST, digits), "(", round(NCT$SE1110_EST, digits), ")", sep="")
    eff4 <- paste(round(NCT$EFF0100_EST, digits), "(", round(NCT$SE0100_EST, digits), ")", sep="") 
  }
  else {
    eff1 <- paste(round(NCT$EFF1000_EST, digits), sep="")
    eff2 <- paste(round(NCT$EFF1101_EST, digits), sep="")
    eff3 <- paste(round(NCT$EFF1110_EST, digits), sep="")
    eff4 <- paste(round(NCT$EFF0100_EST, digits), sep="")  
  }
  
  # Calculate percentages
  total_n <- sum(NCT$NOBS)
  percentages <- sprintf("%.1f%%", (NCT$NOBS / total_n) * 100)
  nobs_with_pct <- paste0(NCT$NOBS, " (", percentages, ")")
  
  V(grafo_tree)$labels <- paste(eff1, eff4, nobs_with_pct, sep="\n")
  
  par(mar = margins)
  
  # Create layout with more horizontal spacing
  layout <- layout_as_tree(grafo_tree)
  # Increase horizontal spacing by spreading x-coordinates much more
  layout[, 1] <- layout[, 1] * 3.5
  
  NCTPLOT <- plot(grafo_tree,
                  layout = layout,
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
                  vertex.label.cex = vertex_label_cex * 0.85,
                  vertex.frame.color = vertex_frame_color,
                  edge.color = edge_color,
                  vertex.frame = 2,
                  edge.label = E(grafo_tree)$label,
                  vertex.shape = vertex_shape,
                  vertex.size = vertex_size * 1.4,
                  vertex.size2 = vertex_size2 * 1.4)
  
  legend('topleft', text.font = 4, text.col = legend_color,
         xjust = 0, adj = 0.2, legend = c(expression(paste(tau, "(1,0;0,0)",sep = "")),
                                          expression(paste(tau,"(0,1;0,0)", sep = "")), "N"),
         title = "Legend", title.col  = legend_title_color)
  
  title(title, cex.main = main_cex, col.main = main_color,
        adj = adj, font.main = main_font)
}