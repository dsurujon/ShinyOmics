library(igraph)

make_node_table<-function(edgetable){
  #make new nodes table from edges
  nodes<-unique(c(edgetable$source,edgetable$target))
  nodeIDs<-c(1:(length(nodes)))
  nodetable<-data.frame("label"=nodes,"id"=nodeIDs)
  
  #add network stats
  mynet <- graph_from_edgelist(as.matrix(edgetable[,c('source','target')]))
  nodetable['Degree'] <- degree(mynet, v=nodetable$label)
  nodetable['Betweenness'] <- betweenness(mynet, v=nodetable$label)
  eigencen <- eigen_centrality(mynet)$vector
  nodetable['Eigencentrality'] <- eigencen[match(nodetable$label, names(eigencen))]
  
  return(nodetable)
}
