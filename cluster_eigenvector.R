var.intersect <- intersect( df.E12.5@var.genes, df.E14.5@var.genes )

temp <- sapply( unique( df0.E12.5@meta.data$res.1 ), function( xx ){
  temp <- df0.E12.5@scale.data[ df0.markers$gene[1:1000],df0.E12.5@meta.data$res.1 == xx ]
  prcomp(temp)$x[ ,1 ]
  })
# temp <- temp[ var.intersect, ]

temp.2 <- sapply( unique( df0.E14.5@meta.data$res.1 ), function( xx ){
  temp <- df0.E14.5@scale.data[ df0.markers$gene[1:1000],df0.E14.5@meta.data$res.1 == xx ]
  prcomp(temp)$x[ ,1 ]
  })
# temp.2 <- temp.2[ var.intersect, ]

res <- apply( temp, 2, function( xx ){apply( temp.2, 2, function( yy ){cor( xx, yy, method = "s" )} )} )
#res[ res < 0.2 ] <- NA
heatmap.3( res, Rowv = T, Colv = T, na.color = "black" )

ddrtree_res <- DDRTree( temp, dimensions = 2 )
mm@reducedDimS <- ddrtree_res$Z
mm@reducedDimK <- ddrtree_res$Y
adjusted_K <- t( reducedDimK( mm ) )
dp <- as.matrix( dist( adjusted_K ) )
mm@cellPairwiseDistances <- dp
gp <- graph.adjacency( dp, mode = "undirected", weighted = TRUE )
dp_mst <- minimum.spanning.tree( gp )
mm@minSpanningTree <- dp_mst
mm@dim_reduce_type <- "DDRTree"

cc_ordering <- extract_ddrtree_ordering( mm, root_cell = NULL )
pData(mm)$Pseudotime <- cc_ordering$pseudo_time
pData(mm)$Pseudodev <- cc_ordering$pseudo_time

ggplot( pData(mm), aes( x = Pseudotime, y = Pseudodev, color = res@meta.data$time_point ) ) +
  geom_point() +
  scale_color_discrete( name = "" )

# mm <- setOrderingFilter(mm, ordering_genes = df@var.genes)

S_matrix <- reducedDimS( mm )
data_df <- data.frame( t( S_matrix[1:2,] ) )
colnames(data_df) <- c( "data_dim_1", "data_dim_2" )
rownames(data_df) <- res@cell.names#df@cell.names
data_df <- merge( data_df, pData( mm ), by.x = "row.names", by.y = "row.names" )

dp_mst <- minSpanningTree( mm )
edge_list <- as.data.frame( get.edgelist( dp_mst ) )
colnames(edge_list) <- c( "source", "target" )

reduced_dim_coords <- reducedDimK( mm )
ica_space_df <- data.frame( Matrix::t( reduced_dim_coords[1:2, ] )    )
colnames(ica_space_df) <- c( "prin_graph_dim_1", "prin_graph_dim_2" )
ica_space_df$sample_name <- row.names( ica_space_df )

edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
edge_df <- plyr::rename( edge_df, c( prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "source_prin_graph_dim_2" ))
edge_df <- merge(edge_df, ica_space_df[, c( "sample_name", "prin_graph_dim_1", "prin_graph_dim_2" )], 
                 by.x = "target", by.y = "sample_name", all = TRUE )
edge_df <- plyr::rename( edge_df, c( prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "target_prin_graph_dim_2" ) )
