#' @title
#' Plot of a Network Causal Tree

#' @description
#' Plots the Network Causal Tree results

#' General
#' @param NCT NetworkCausalTrees results.
#' @param varnames (Ordered) covariates names vector.
#' @param output If "detection" only point estimates are reported, if
#' "estimation" both estimated effects and variances are reported.
#' @param ot TRUE if the results. includes an optimal treatment evaluation,
#' FALSE otherwise (default: FALSE).
#' Edges
#' @param ewidth Edge width.
#' @param elabelcex Edge's label cex.
#' @param edge.arrow.size Edge's arrow size.
#' @param ecolor Edges color.
#' @param elabelfamily Edge label family.
#' @param elabelcolor Edges' label color.
#' Vertices
#' @param vsize Vertex size.
#' @param vsize2 Vertex size (for rectangular vertices).
#' @param vshape Vertex shape.
#' @param vlabelfont Vertex label font.
#' @param vlabelcex Vertex Label cex.
#' @param coloreff Estimated Effect according to which nodes are colored
#' It can be "1000", "0100","1110","1101".
#' @param vcolor String Vector of four colors: the first two identify the
#' extremes of the color palette employed to depict negative values, while the
#' second two identify the extremes of the color palette employed to depict
#' positive values. Note that the intensity of the color is related to the
#' Statistic of the estimated effect in each leaf (Default values represent
#' negative estimates with greens and positive estimates with blues).
#' @param vframecolor Vertex's frame color.
#' @param vlabelcolor Vertex label color.
#' Title
#' @param title Title.
#' @param cex.main Title cex.
#' @param adj adj.
#' @param col.main Title color.
#' Legend
#' @param colleg Legend color.
#' @param coltitleleg Legend's title color.
#'
#' @return ADD HERE A PLOT DESCRIPTION
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
                     vcolor = c("seagreen4","seagreen1","lightblue1","dodgerblue2"),
                     vsize = 32,
                     vsize2 = 32,
                     vshape = "circle",
                     vframecolor = "black",
                     vlabelcolor = "black",
                     vlabelcex = 0.8,
                     vlabelfont = 1,
                     elabelcex = 0.7,
                     elabelfamily = "sans",
                     ecolor = "black",
                     ewidth = 0.3,
                     edge.arrow.size = 0.5,
                     elabelcolor = "black",
                     ot = FALSE,
                     colleg = "black",
                     coltitleleg = "black",
                     font.main = 1,
                     cex.main = 1,
                     adj = 1,
                     col.main = "black"){

  options(warn=-1)
  NCT$NOBS<-NCT$NOBS_EST+NCT$NOBS_TR

  for (p in 1:length(varnames)) {
    NCT$FILTER<-gsub(pattern = paste("X.",p,sep=""),replacement = varnames[p],x=as.character(NCT$FILTER) )
  }

  NCT$FILTER<-gsub(pattern = "data_tree",replacement ="",x=as.character(NCT$FILTER) )
  NCT$FILTER<-gsub(pattern = "[$]",replacement ="",x=as.character( NCT$FILTER ) )
  NCT$FILTER<-gsub(pattern = "_bin",replacement ="",x=as.character( NCT$FILTER ) )
  NCT$FILTER[which(NCT$FILTER!="NA")]<-paste0(" & ", NCT$FILTER[which(NCT$FILTER!="NA")])
  NCT$FILTER[which(NCT$FILTER!="NA")]<-paste0(" NA ",NCT$FILTER[which(NCT$FILTER!="NA")])

  tree_data<-as.Node(NCT, mode = "table",
                     pathName = "FILTER", pathDelimiter = " & ", colLevels = NULL,
                     na.rm = TRUE)
  lati_tree<-ToDataFrameNetwork(tree_data)
  lati_tree<-as.matrix(lati_tree)
  nomelati<-NULL
  for(i in 1:nrow(lati_tree)){
    nomelati_i<-tail(strsplit(lati_tree[,2], "/")[[i]],n=1)
    nomelati_i<-gsub(pattern = "X.",replacement="",x=nomelati_i)
    nomelati=c(nomelati,nomelati_i)
  }
  grafo_tree<-graph_from_edgelist(lati_tree, directed = TRUE)
  V(grafo_tree)$TAU1000<-NCT$EFF1000_EST
  V(grafo_tree)$TAU1101<-NCT$EFF1101_EST
  V(grafo_tree)$TAU1110<-NCT$EFF1110_EST
  V(grafo_tree)$TAU0100<-NCT$EFF0100_EST

  if(output=="estimation"){
    V(grafo_tree)$SE1000<-NCT$SE1000_EST
    V(grafo_tree)$SE1101<-NCT$SE1101_EST
    V(grafo_tree)$SE1110<-NCT$SE1110_EST
    V(grafo_tree)$SE0100<-NCT$SE0100_EST

    if(coloreff=="1000"){
      V(grafo_tree)$stat<- V(grafo_tree)$TAU1000/ V(grafo_tree)$SE1000
    }
    if(coloreff=="1101"){
      V(grafo_tree)$stat<- V(grafo_tree)$TAU1101/ V(grafo_tree)$SE1101
    }
    if(coloreff=="0100"){
      V(grafo_tree)$stat<- V(grafo_tree)$TAU0100/ V(grafo_tree)$SE0100
    }
    if(coloreff=="1110"){
      V(grafo_tree)$stat<- V(grafo_tree)$TAU1110/ V(grafo_tree)$SE1110
    }

    cols1 = colorRampPalette(c(vcolor[1],vcolor[2]))(4)
    cols2 = colorRampPalette(c(vcolor[3], vcolor[4]))(4)
    cols<-c(cols1,cols2)

    V(grafo_tree)$color<-cols[1]
    V(grafo_tree)$color[which(V(grafo_tree)$stat< (-2.58))]<-rep(cols[1], length(which(V(grafo_tree)$stat< (-2.58))))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>= (-2.58) & V(grafo_tree)$stat<(-1.96) )]<-rep(cols[2], length(which(V(grafo_tree)$stat>= (-2.58) & V(grafo_tree)$stat<(-1.96) )))
    V(grafo_tree)$color[which(V(grafo_tree)$stat >= (-1.96) & V(grafo_tree)$stat<(-1.66) )]<-rep(cols[3], length(which(V(grafo_tree)$stat>= (-1.96) & V(grafo_tree)$stat<(-1.66) )))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>= (-1.66) & V(grafo_tree)$stat<(0 )) ]<-rep(cols[4], length(which(V(grafo_tree)$stat>= (-1.66) & V(grafo_tree)$stat<(0 ))))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>= 0 & V(grafo_tree)$stat<(1.66 ))]<-rep(cols[5],length(which(V(grafo_tree)$stat>= 0 & V(grafo_tree)$stat<(1.66 ))))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>=1.66 & V(grafo_tree)$stat<1.96 )]<-rep(cols[6],length(which(V(grafo_tree)$stat>=1.66 & V(grafo_tree)$stat<1.96 )))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>=1.96 & V(grafo_tree)$stat<2.58 )]<-rep(cols[7],length(which(V(grafo_tree)$stat>=1.96 & V(grafo_tree)$stat<2.58 )))
    V(grafo_tree)$color[which(V(grafo_tree)$stat>=2.58) ]<-rep(cols[8],length(which(V(grafo_tree)$stat>=2.58)))
  }

  if(output=="detection"){
    V(grafo_tree)$color<-vcolor[1]
    print("Unique color applied because NCT does not provide standard errors")
  }

  latinomefine<-rep(NULL,length(nomelati))
  edgelabel<-vector(mode = "list", length = length(nomelati))
  for (l in 1:length(nomelati)){
    if(length(grep(x=nomelati[l],pattern=">="))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[>=])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[>=])",perl = TRUE)[[1]][2:(lengths(gregexpr(",", strsplit(nomelati[l],"(?=[>=])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern="<="))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[<=])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[<=])",perl = TRUE)[[1]][2:(lengths(gregexpr(",", strsplit(nomelati[l],"(?=[<=])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern="<="))==0 & length(grep(x=nomelati[l],pattern="<"))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[<])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[<])",perl = TRUE)[[1]][2:(lengths(gregexpr(",", strsplit(nomelati[l],"(?=[<])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern=">="))==0 & length(grep(x=nomelati[l],pattern=">"))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[>])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[>])",perl = TRUE)[[1]][2:(lengths(gregexpr(",", strsplit(nomelati[l],"(?=[>])",perl = TRUE))) + 1)],collapse = " ")
    }
    latinomefine[l]<-paste(edgelabel[[l]][1],gsub("\\s", "", x=edgelabel[[l]][2]),sep="\n")

  }
  latinomefine<-gsub(pattern = "<\\(1\\)",replacement ="=0",x=as.character(latinomefine) )
  latinomefine<-gsub(pattern = ">=\\(1\\)",replacement ="=1",x=as.character(latinomefine) )

  E(grafo_tree)$label<-latinomefine

  if (output=="estimation") {
    eff1<-paste(round(NCT$EFF1000_EST,2),"(",round(NCT$SE1000_EST,2),")",sep="")
    eff2<-paste(round(NCT$EFF1101_EST,2),"(",round(NCT$SE1101_EST,2),")",sep="")
    eff3<-paste(round(NCT$EFF1110_EST,2),"(",round(NCT$SE1110_EST,2),")",sep="")
    eff4<-paste(round(NCT$EFF0100_EST,2),"(",round(NCT$SE0100_EST,2),")",sep="")
  }
  if (output=="detection") {
    eff1<-paste(round(NCT$EFF1000_EST,2),sep="")
    eff2<-paste(round(NCT$EFF1101_EST,2),sep="")
    eff3<-paste(round(NCT$EFF1110_EST,2),sep="")
    eff4<-paste(round(NCT$EFF0100_EST,2),sep="")
  }
  if (ot==TRUE) {
    ot<-paste("(", as.numeric(round(NCT$OT,2)),")",sep="") #";",as.numeric(round(NCT$OT_NODES_BC,2)),")",sep = "")
    V(grafo_tree)$labels<-paste(eff1,eff4,NCT$NOBS,ot,sep="\n")
  } else {
    V(grafo_tree)$labels<-paste(eff1,eff4,NCT$NOBS, sep="\n")
  }

  NCTPLOT<-plot(grafo_tree,layout=layout_as_tree(grafo_tree), edge.label.color=elabelcolor,
                edge.width=ewidth,edge.label.cex=elabelcex, edge.label.family=elabelfamily, vertex.color=V(grafo_tree)$color,
                vertex.label.dist = 0, vertex.label.color=vlabelcolor,
                vertex.label=V(grafo_tree)$labels,  vertex.label.font=vlabelfont, vertex.label.cex=vlabelcex,
                vertex.frame.color=vframecolor, edge.color=ecolor, vertex.frame=2, edge.arrow.size=edge.arrow.size,
                edge.label=E(grafo_tree)$label,vertex.shape=vshape,vertex.size=vsize,vertex.size2=vsize2)

  if (ot==TRUE) {
    legend('topleft', text.font = 4, text.col=colleg, xjust=0, adj=0.2,legend=c(expression(paste(tau,"(1,0;0,0)",sep="")),
                                                                                expression(paste(tau,"(0,1;0,0)",sep="")),"N","OT"),
           title = "Legend",title.col =coltitleleg)
  } else {
    legend('topleft', text.font=4,text.col=colleg, xjust=0, adj=0.2,legend=c(expression(paste(tau,"(1,0;0,0)",sep="")),
                                                                             expression(paste(tau,"(0,1;0,0)",sep="")),"N"),
           title = "Legend",title.col =coltitleleg)
  }
  if(output=="estimation"){
  legend('topright', legend=c("z < -2.58", " -2.58 <= z < -1.96  ",
                              " -1.96 <= z < -1.66  ", " -1.66 <= z < 0  ",
                              " 0 <= z < 1.66 ", "1.66 <= z < 1.96", "1.96 <= z < 2.58","2.58 < z" ),
         col=cols, pch=15, cex=0.7,
         text.font=4, bg='white')}

  title(title,cex.main=cex.main,col.main=col.main,adj=adj,font.main=font.main)
  invisible(NCTPLOT)
}

