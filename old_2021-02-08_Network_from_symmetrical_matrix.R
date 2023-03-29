# Update 2021-02-05: v1.2 I added a new variant of the clustering function in case the input graph is too large (results in memory bloating with segment dump). This new version doesn't carry out edge_betweenness calculation (both the calculation per se and the clustering approach) and the spinglass algorithm, as these are the most memory-intensive
# Update 2021-02-04: v1.1 I added a new input parameter: print_names (boolean), to detect if the user wants to remove labels (if FALSE). In this case, a cross reference table is printed at the beggining for traceability, then nodes are renamed as v1, v2... vn
# Started on 2021-01-18 (v1)
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock" 
# It is intended to carry out a network analysis based on input adjacency matrix (symmetrical)
# The input matrix is expected to be passed as stdin so that it can be used as part of a pipeline
# The first row and column should match, containing feature names in the same order
# A symmetrical matrix is expected (e.g. distance or correlation matrix), using a linear scale.
# Both, positive and negative values are allowed but only matrices designated as "correlations" are considered for dual-colored edges. Others having both positive and negative values are shifted so that all items are >= 0, as centrality calculations requiere positive-only weights (no further transformations beyond this recenter step are carried out with the raw data)
# Network construction is already based on an adjacency matrix but can carry an optional (the simplest) sparsification step consisting on filtering based on a user-defined cutoff minimum value (in case correlation values are selecte, this will be sign-aware)
# Networks are build using edge maps that may be filtered through sparsification, so, removing all edges connecting a node will result on the node being deleted from the resulting network
# IMPORTANT NOTE 1: This script asumes input data reflect closeness (the further from 0, the strongest the relationship between nodes). Distances and semi-metrics from vegan (like the ones QIIME uses) are calculated as similarities, meaning 1 is the max similarity whereas 0 is totally unrelated (similarly to sign-less correlation values). The script is intended for this type of differences, and should be adjusted accordingly if weights mean refer to euclidean or similar distances (where larger weights are more distant nodes).
# IMPORTANT NOTE 2: Weights are considered linear for edge weight as they will be scaled to a range of 0-1 for some calculations. Thus, if log values are present, please convert them to linear scale before using this script.
# Max limit: The script should hold for at least 700 vertices and 3500 edges (it hasn't been tested for more).

# Run as follows:
# cat symmetrical_matrix.tsv|Rscript Network_from_symmetrical_matrix.R <metadata_input> <prefix_output> <gt0_if_correlation_matrix> [min_value (optional)]
# Tested with command:
# cat /home/rod/Documents/02_Collaborations/Geneticas/01_input_tables/weighted_unifrac_dm.tsv|Rscript Network_from_symmetrical_matrix.R /home/rod/Documents/02_Collaborations/Geneticas/01_input_tables/metadata.tsv test FALSE 0.5
# Internally: 
#  df <- as.matrix(read.table("01_input_tables/weighted_unifrac_dm.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F))
#  meta <- "/home/rod/Documents/02_Collaborations/Geneticas/01_input_tables/metadata.tsv"
#  print_names <- FALSE
#  is_corr <- as.logical(0)
#  prefix <- "test"
#  min_val <- 0
# df2 <- df*matrix(sample(c(-1,1),35*35,replace=T),35,35) # simulate correlation matrix (-1 to 1)
# df <- as.matrix(read.table("03_corr_tables/genetics_lvl2-H_grp_A-corr_pearson_rows.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F))
# meta <- "01_input_tables/meta_all.tsv"
# is_corr <- TRUE
# min_val <- 0.5
# prefix <- "test2"
# print_names <- TRUE


### PRE-LOAD FILES AND PARAMETERS ###
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, two arguments should be included: <metadata_table.tsv> <prefix_output> <TRUE/FALSE_if_correlation_matrix> <TRUE/FALSE_to_print_names> [min_value (optional)]
  stop("A minimum of 4 arguments are mandatory: cat input_matrix.tsv|Rscript Network_from_symmetrical_matrix.R <metadata_table.tsv> <prefix_output> <TRUE/FALSE_if_correlation_matrix> <TRUE/FALSE_to_print_names> [min_value (optional)]", call.=FALSE)
}
meta <- as.character(args[1])  # Get a path and file for the metadata table (tsv)
prefix <- as.character(args[2]) # Get a string handle to create output names
is_corr <- as.logical(args[3]) # Boolean flag for indictaing the input is a correlation matrix (-1 to 1)
print_names <- as.logical(args[4]) # Boolean flag indicating whether to print full vertex names or not
min_val <- as.numeric(args[5]) # Optional, get the minimum value cutoff. Adjust this according to the value range in the matrix. This ignores signs if is_corr is passed as TRUE

# Input matrix
df <- as.matrix(read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)) # Load the input matrix from stdin
if(isSymmetric(df)==FALSE){stop(paste0("Aborting: Matrix is not symmetrical."), call.=FALSE)} #  abort if the input matrix does not comply with being symmetrical
if(!file.exists(meta)){stop(paste0("Aborting: Metadata file not found. Please check if file name and path are correct."), call.=FALSE)} # Abort if the metadata file is not present
# Input metadata
mdf <- data.frame(read.table(meta, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F, stringsAsFactors=TRUE)) # Load the metadata companion matrix. All nodes must be present
if(sum(is.na(match(rownames(df),rownames(mdf))))>0) {stop(paste0("Aborting: Samples are missing from the metadata file. Please check if all samples are present."), call.=FALSE)} # Abort if the metadata file is not present
temp <- colnames(mdf)
mdf <- as.data.frame(mdf[rownames(df),]) # subset the metadata to match the input matrix order
rownames(mdf) <- rownames(df); colnames(mdf) <- temp # This is a fix in case only one column is present as the last step will try to treat it as a factor vector
if(print_names==FALSE){ # If the user requested not to use the original names
	new_names <- rownames(df)
	names(new_names) <- paste0("v",1:nrow(df))
	write.table(as.data.frame(new_names),paste0(prefix,"-cut-",min_val,"-vertices_names.tsv"), sep="\t", quote=FALSE, col.names=FALSE)
	rownames(df) <- names(new_names) # After saving the original names for cross reference, rename all vertices
	colnames(df) <- names(new_names)
	rownames(mdf) <- names(new_names)
}
if(is.na(min_val)){min_val=min(df[unlist(df)>0])} # Default value for min_val (an optional argument) is set to the smallest value in the table (i.e. no filters)
print("Loaded parameters:")
print(paste("Metadata file name:",meta))
print(paste("Prefix:",prefix))
print(paste("Correlation:",is_corr))
print(paste("Minimum value cutoff:",min_val))

### LOAD LIBRARIES AND FUNCTIONS ###
# Libraries
library("igraph")
library("reshape2")
library("RColorBrewer") # display.brewer.all() # these are available colors schemes in 

#Functions
rescale_vect <- function(vect,init,end){
	if(length(table(vect))==1){
		out <- rep(init,length(vect))
	}else{
		out <- (end-init)/(max(vect)-min(vect))*(vect-min(vect))+init
	}
	return(out)
}
load_edge_map <- function(sym_mat){
	sym_mat[lower.tri(sym_mat, diag=T)] <- NA # remove the lower triangle of the symmetrical matrix
	edge_map <- na.omit(melt(sym_mat, value.name ="weight")) # get long format for all edges. The lower triangle and the diagonal are ignored
	row.names(edge_map) <- 1:nrow(edge_map) # reset edge numbers
	return(edge_map)
}
filter_map <- function(edge_map,corr,cutoff){
	edge_map <- cbind(edge_map, "org_weight"=edge_map[,3]) # First, make a copy of the original weights
	if(corr){ # If it is a correlation matrix (this flag is specified by the user)
		edge_map <- edge_map[abs(edge_map[,3])>=cutoff,] # Then consider the absolute value when filtering
		edge_map[,3] <- abs(edge_map[,3]) # and change weights to positive values
	}else{ # if central value is not 0 (flagged as correlations)
		edge_map <- edge_map[edge_map[,3]>=cutoff,] # In any other case, just cut any values below the preset cutoff
		if(min(edge_map[,3])<0){ # now, if there are still any numbers below 0
			edge_map[,3] <- edge_map[,3]+abs(min(edge_map[,3]))# we need to shift values to start at 0 (sum the absolute of the min value to make them 0)
		}
	}
	edge_map <- cbind(edge_map, "scaled_01"=rescale_vect(edge_map[,3],0,1)) # now scale the weight vector to 0-1 and append as a column (it uses a custom function, defined above)
	return(edge_map)
}
set_node_col <- function(n){ # set a color scheme based on how many factors may be present
	if(n <= 2){ # Min allowed for the color scheme is 2, this adds an exception for less
		out <- brewer.pal(3,"Set2")[1:n]
	} else if(n <= 8) { # default scheme is used if less than 9 (max is 8 for "Set2" color scheme)
		out <- brewer.pal(n,"Set2")
	} else { # if more are required, a gradient is used instead, based on the "Set2" scheme
		out <- colorRampPalette(brewer.pal(8,"Set2"))(n)
	}
	return(out)
}
prune_net <- function(in_net){
	deg <- degree(in_net, mode="all") # get not degrees (since it is undirected, we use both in-degree and out-degree values)
	print(paste0("Found ",sum(degree(in_net)>0), " connected nodes"))
	print(paste0("Ignoring ",sum(degree(in_net)==0), " unconnected nodes"))
	lone <- which(degree(in_net)==0) # detect disconnected nodes and remove them
	in_net = delete.vertices(in_net, lone)
	return(in_net)
}

structure_analysis <- function(in_net, pre){
	meta <- mdf[V(in_net),] # add the metadata for plotting colors (subset only the existing nodes
	nodes <- gorder(in_net)
	edges <- gsize(in_net)
	filt_mean_weight <- mean(E(in_net)$org_weight)
	filt_weight_sd <- sd(E(in_net)$org_weight)
	filt_weight_median <- median(E(in_net)$org_weight)
# 	quant <- quantile(E(in_net)$org_weight,seq(0,1,0.1)) # too long for an output
	dens <- edge_density(in_net, loops=F) # Get the proportion of present edges from all possible edges (0-1)
	deg <- degree(in_net, mode="all") # get not degrees (since it is undirected, we use both in-degree and out-degree values)
	mean_dist <- mean_distance(net, directed=F)
# 	dist_mat <- distances(in_net) # Optionally, get the total distance for each two nodes (redundant)
	cliques <- largest_cliques(in_net) # get the largest clique (totally connected subgraph)
	cliques <- unlist(lapply(cliques,function(x) paste(V(in_net)$name[x],collapse='<->'))) # create a simple vector
	names(cliques) <- paste("Max clique",1:length(cliques))
	trans_g <- transitivity(in_net, type="globalundirected") # Calculate the clustering/transitivity coeficient. This is the ratio of triangles (triads) to edges (pairs) in the network (all permutations
	trans_l <- transitivity(in_net, type="localundirected") # same for every node
	diam <- diameter(in_net) # Get the max geodesic distance in the whole in_network (considering the shortest between nodes)
	diam_nodes <- get_diameter(in_net) # Get the corresponding nodes that have this diameter
	# E(in_net, path=diam_nodes) # This lists the actual paths
	vertcol <- rep("cornflowerblue", vcount(in_net))
	vertcol[diam_nodes] <- "coral1"
	edgecol <- rep("gray80", ecount(in_net))
	edgecol[E(in_net, path=diam_nodes)] <- "orange" 
	deg_dist <- degree_distribution(in_net, cumulative=T, mode="all")
	# Calculate a centrality matrix (multiple centrality measurements)
	eigen <- centr_eigen(in_net, directed=T, normalized=F) # only this one needs to be filtered
	tempnet <- in_net # This is a fix introduced as betweenness and closeness is calculated as actual node distances. By default, these two centrality mentions are unweighted.
	E(tempnet)$weight <- 1.001-E(tempnet)$scaled_01 # adjust so that they are all positive distances (they are not the original values but proportional distances are conserved)
	centrality <- round(cbind("degree_u"=unlist(centr_degree(in_net, mode="all", normalize=F)),"eigenvector_u"=c(eigen$vector,eigen$centralization,eigen$theoretical_max),"closeness_u"=unlist(centr_clo(in_net, mode="all", normalize=F)),"betweenness_u"=unlist(centr_betw(tempnet, directed=FALSE, normalized=F))),2) # Calculate other centrality measurements and create a dataframe
	centrality <- cbind(centrality, "transitivity_u"=c(trans_l, NA, trans_g)) # append the previously calculated transitivity
	centrality <- cbind(centrality, "closeness_w" = c(closeness(tempnet, weights=1.001-E(tempnet)$scaled_01),NA,NA))
	centrality <- cbind(centrality, "betweenness_w" = c(betweenness(tempnet, weights=1.001-E(tempnet)$scaled_01),NA,NA))
	centrality <- cbind(centrality, "hub.score"=c(hub_score(in_net, weights=NULL)$vector,NA,NA))
# 	centrality <- cbind(centrality, "authority.score"=c(authority_score(in_net, weights=NA)$vector,NA,NA)) # same as hubs if undirected
	centrality <- cbind(centrality, "vertex.strength"=c(graph.strength(in_net),NA,NA))
	rownames(centrality)[1:vcount(in_net)] <- names(deg) # rename the vertex
	# Now for global specs
	global_specs <- data.frame(min_val, nodes, edges, filt_mean_weight, filt_weight_sd, filt_weight_median, dens, diam, mean_dist, "nodes_in_diam"=paste(V(in_net)$name[get_diameter(in_net)], collapse='-'))
	global_specs <- t(cbind(global_specs, t(cliques)))
	write.table(global_specs,paste0(prefix,"-cut-",min_val,"-",pre,"-glob_specs.tsv"), sep="\t", quote=FALSE, col.names=FALSE)
	pdf(paste0(prefix,"-cut-",min_val,"-",pre,"-net.pdf"))
		par(las=1)
		plot(in_net, vertex.color=vertcol, vertex.frame.color="white", vertex.label.color="black", edge.color=edgecol, vertex.size=rescale_vect(deg,15,25), vertex.label.family="Helvetica", main="Degree and diameter (longest path)", layout=lo)
		mtext(paste0("Cutoff value: ", min_val))
		par(mfrow = c(1,2))
		plot(net, vertex.size=rescale_vect(centrality[1:(nrow(centrality)-2),8]*50,15,25), main="Main Hubs (nodes)", vertex.color=vertcol, vertex.frame.color="white", vertex.label.color="black", edge.color=edgecol, vertex.label.family="Helvetica", vertex.label.cex=0.3, layout=lo)#, vertex.label=NA) # Plot the same network with larger hubs (adjust vertex.label.size or add vertex.label=NA accordingly)
		mtext(paste0("Cutoff value: ", min_val))
		plot(net, vertex.size=rescale_vect(centrality[1:(nrow(centrality)-2),9],15,25), main="Vertex strength", vertex.color=vertcol, vertex.frame.color="white", vertex.label.color="black", edge.color=edgecol, vertex.label.family="Helvetica", vertex.label.cex=0.3, layout=lo) # Plot the network with larger authorities
		mtext(paste0("Cutoff value: ", min_val))
		par(mfrow = c(1,1))
		hist(deg, breaks=1:vcount(in_net)-1, main="Node degree",col="cornflowerblue", xlab="Degree of connectivity (undirected)", ylab="Nodes", border="white")
		mtext(paste0("Cutoff value: ", min_val))
		plot(x=0:max(deg), y=1-deg_dist, type='o', pch=19, cex=1.2, col="coral1", xlab="Degree", ylab="Cumulative Frequency", main="Degree distribution")
		mtext(paste0("Cutoff value: ", min_val))
		hist(trans_l, main="Cluster/transitivity coefficient",col="cornflowerblue", xlab="Cluster coefficient (node transitivity)", ylab="Nodes", border="white")
		mtext(paste0("Cutoff value: ", min_val))
		plot(cumsum(sort(trans_l[!is.na(trans_l)])),seq(0,1,1/sum(!is.na(trans_l)))[-1], type='o', pch=19, cex=1.2, col="coral1", xlab="Cluster coefficient (transitivity)", ylab="Cumulative Frequency", main="Cluster/transitivity coefficient distribution")
		mtext(paste0("Cutoff value: ", min_val))
# 		plot(cumsum(trans_l[!is.na(trans_l)]),seq(0,1,1/sum(!is.na(trans_l)))[-1], type='o', pch=19, cex=1.2, col="coral1", xlab="Cluster coefficient (transitivity)", ylab="Cumulative Frequency", main="Cluster/transitivity coefficient distribution")
		par(mfrow = c(2,2)) # Plot histograms of centrality measurments (in a single plot)
		hist(centrality[1:(nrow(centrality)-2),1], main="Node Degree Centrality",col="cornflowerblue", xlab="Node connections (undirected)", ylab="Nodes", border="white")
		hist(centrality[1:(nrow(centrality)-2),2], main="Node Eigenvector Centrality",col="cornflowerblue", xlab="Neighbouring connectivity", ylab="Nodes", border="white")
		hist(centrality[1:(nrow(centrality)-2),3], main="Node Unweighted Closeness Centrality",col="cornflowerblue", xlab="Node closeness", ylab="Nodes", border="white")
		hist(centrality[1:(nrow(centrality)-2),4], main="Node Unweigted Betweenness Centrality",col="cornflowerblue", xlab="Node betweenness", ylab="Nodes", border="white")
		hist(centrality[1:(nrow(centrality)-2),5], main="Node Weighted Closeness Centrality",col="cornflowerblue", xlab="Node closeness", ylab="Nodes", border="white")
		hist(centrality[1:(nrow(centrality)-2),6], main="Node Weigted Betweenness Centrality",col="cornflowerblue", xlab="Node betweenness", ylab="Nodes", border="white")
		par(mfrow=c(1,1))
	dev.off()
	if(print_names==FALSE){ # If the user requested not to use the original names, add the names before print
		centrality <- as.data.frame(cbind(centrality, "vertex_name"=new_names[rownames(centrality)])) # this uses a simple trick to keep the names of df "centrality"
		centrality <- cbind(centrality[ncol(centrality)],centrality[-ncol(centrality)]) # now shuffle columns into the correct order
	}
	write.table(centrality,paste0(prefix,"-cut-",min_val,"-",pre,"-node_specs.tsv"), sep="\t", quote=FALSE, col.names=NA)
	return(list(centrality, global_specs))
}
# structure_analysis(net,"full")
plot_meta <- function(in_net, metadata, pre){
	range <- round(seq(-1,1,0.01),2) # Create a range for expected rho values (201 items, from -1 to 1 by default)
	temp <- colnames(metadata)
	metadata <- as.data.frame(metadata[V(in_net)$name,]) # subset the metadata table in case some nodes were unconnected
	rownames(metadata) <- V(in_net)$name; colnames(metadata) <- temp # This fix deals with single-column metadata files
	dict <- colorRampPalette(c("red", "white", "blue"))(length(range)) # and create a dictionary for colors
	names(dict) <- range
	if(is_corr){ # if the user specified weights have values between -1 and 1 (correlation-like)
		edge_cols <- dict[as.character(round(E(in_net)$weight*E(in_net)$org_weight/abs(E(in_net)$org_weight),2))]
	}else{ # in any other case, use values standardized between 0 and 1
		edge_cols <- dict[as.character(round(E(in_net)$scaled_01,2))]
	}
# 	lo <- layout_with_fr(in_net) # Create a fixed layout for all Fruchterman and Reingold graphs
	pdf(paste0(prefix,"-cut-",min_val,"-",pre,"-groups.pdf"))
	for(i in 1:ncol(metadata)){
		node_cols <- set_node_col(1) # By default, set to 1 color (if no factors are present)
		sizes <- 15 # set default size to 25
		items <- metadata[,i] # Load each vector
		temp <- rep(1,length(items)) # By default, use the first color
		if(is.factor(items)){ # if it is a factor, create a new color scheme
			node_cols <- set_node_col(length(levels(items)))
			temp <- as.numeric(items)
		}else if(is.numeric(items)){
			sizes <- rescale_vect(items,10,25)
		}
		plot(net, vertex.color=sapply(temp, function(x) node_cols[x]), vertex.frame.color="white", vertex.label.color="black", vertex.label.family="Helvetica", main=paste0("Network depicting group: ", names(metadata)[i]), edge.width=rescale_vect(E(in_net)$weight,0.5,3), vertex.size=sizes, edge.color=edge_cols, layout=lo)
		mtext(paste0("Cutoff value: ", min_val))
		plot(net, vertex.color=sapply(temp, function(x) node_cols[x]), vertex.frame.color="white", vertex.label.color="black", vertex.label.family="Helvetica", main=paste0("Network depicting group: ", names(metadata)[i]), edge.width=rescale_vect(E(in_net)$weight,0.5,3), vertex.size=sizes, edge.color=edge_cols,layout=layout.circle(in_net, order=order(items)))
		mtext(paste0("Cutoff value: ", min_val))
	}
	dev.off()
}
# plot_meta(net, mdf, "full")
community_detection <- function(in_net, pre){ # This creates clusters with different algorithms large networks should avoid doing all of them, use community_detection_for_large_graph instead.
# 	lo <- layout_with_fr(in_net) # uncomment if a new plotting scheme is required
	E(in_net)$edge_betw_u <- edge_betweenness(in_net, directed=FALSE, weights=NULL)
	E(in_net)$edge_betw_w <- edge_betweenness(in_net, directed=FALSE, weights=E(in_net)$scaled_01+0.001)
	write.table(as_long_data_frame(in_net),paste0(prefix,"-cut-",min_val,"-",pre,"-edge_specs.tsv"), sep="\t", quote=FALSE, col.names=NA) # Print the edge table (comment this if required as this way too large for big networks)
	deg <- degree(in_net) # recalculate degree
	clu_louv <- cluster_louvain(in_net) # cluster with louvain algorithm requires only positive values for edge weights
# 	clu_edg_btw <- cluster_edge_betweenness(in_net, directed=F, weights=1.01-E(in_net)$scaled_01) # cluster based on edge betweenness (Newman-Girvan).
	# IMPORTANT NOTE: In weighted graphs, this produces a warning, since weights are considered distances for edge betweenness calculations (larger distances mean far apart nodes). According to the warning, modularity calculations with weighted might not make sense since they treat them instead as similarities. This was confirmed through tests so the calculations (using simmilarities) are carried out separately from the actual graphs (using distances) (2021-01-01). To this date, there is no easy fix and using weights is better avoided
	clu_edg_btw <- cluster_edge_betweenness(in_net, directed=FALSE, weights=NULL) # NOTE: This is currently the only way to calculate cluster_edge_betweenness without fault calculations due to the reasons explained above
	clu_lab_prop <- cluster_label_prop(in_net) # clustering based on propagating labels (most common edge paths)
	clu_fast_greed <- cluster_fast_greedy(as.undirected(in_net)) # greedy clustering alogorithm
	clu_walktrap <- walktrap.community(in_net,steps=5) # the original article recommends 5 steps
	clu_spinglass <- spinglass.community(in_net)
	clu_eigen <- leading.eigenvector.community(in_net)
	# Append edge betweennes
	clusters <- as.data.frame(cbind("louvain"=membership(clu_louv), "spinglass"=membership(clu_spinglass), "lab_prop"=membership(clu_lab_prop), "edge_btw_u"=membership(clu_edg_btw), "greedy"=membership(clu_fast_greed), "walktrap"=membership(clu_walktrap), "eigen"=membership(clu_eigen)))
	# Create clusters table
	clusters[nrow(clusters)+1,] <- c(modularity(clu_louv), modularity(clu_spinglass), modularity(clu_lab_prop), modularity(clu_edg_btw), modularity(clu_fast_greed), modularity(clu_walktrap), modularity(clu_eigen))
	rownames(clusters)[nrow(clusters)] <- "modularity"
	if(print_names==FALSE){ # If the user requested not to use the original names
		clusters <- as.data.frame(cbind(clusters,"vertex_name"=new_names[rownames(clusters)])) # this uses a simple trick to keep the names of df "clusters"
		clusters <- cbind(clusters[ncol(clusters)],clusters[-ncol(clusters)]) # now shuffle columns into the correct order
	}
	write.table(clusters,paste0(prefix,"-cut-",min_val,"-",pre,"-clusters.tsv"), sep="\t", quote=FALSE, col.names=NA) # 
	pdf(paste0(prefix,"-cut-",min_val,"-",pre,"-clust.pdf"))
		plot(clu_louv, in_net, vertex.size=sqrt(deg),edge.width=E(in_net)$weight, main="Network clusters based on louvain algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)#, edge.width=rescale_vect(E(in_net)$weight,0.5,3)) #, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_louv$membership)), " / Modularity score: ", round(modularity(clu_louv),4)))
		plot(clu_spinglass, in_net, vertex.size=sqrt(deg), main="Network clustering based on spinglass algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_spinglass$membership)), " / Modularity score: ", round(modularity(clu_spinglass),4)))
		plot(clu_lab_prop, in_net, vertex.size=sqrt(deg), main="Network clustering based on propagating labels algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)#, edge.width=rescale_vect(E(in_net)$weight,0.5,3))#, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_lab_prop$membership)), " / Modularity score: ", round(modularity(clu_lab_prop),4)))
		plot(clu_edg_btw, in_net, vertex.size=sqrt(deg), main="Network clustering based on edge betweenness algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)#, edge.width=rescale_vect(E(in_net)$weight,0.5,3))#, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_edg_btw$membership)), " / Modularity score: ", round(modularity(clu_edg_btw),4)))
		if((length(table(clu_edg_btw$membership))>1) & (length(table(clu_edg_btw$membership))/clu_edg_btw$vcount!=1)){dendPlot(clu_edg_btw, mode="hclust")} #This avoids the error of trying to create a single community graph
		plot(clu_fast_greed, as.undirected(in_net), vertex.size=sqrt(deg), main="Network clustering based on fast greedy algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo) # , edge.width=rescale_vect(E(in_net)$weight,0.5,3))#, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_fast_greed$membership)), " / Modularity score: ", round(modularity(clu_fast_greed),4)))
		if((length(table(clu_fast_greed$membership))>1) & (length(table(clu_fast_greed$membership))/clu_fast_greed$vcount!=1)){dendPlot(clu_fast_greed, mode="hclust")}
		plot(clu_walktrap, in_net, vertex.size=sqrt(deg), main="Network clustering based on walktrap algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_walktrap$membership)), " / Modularity score: ", round(modularity(clu_walktrap),4)))
		if((length(table(clu_walktrap$membership))>1) & (length(table(clu_walktrap$membership))/clu_walktrap$vcount!=1)){dendPlot(clu_walktrap, mode="hclust")}
		plot(clu_eigen, in_net, vertex.size=sqrt(deg), main="Network clustering based on leading eigenvector algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_eigen$membership)), " / Modularity score: ", round(modularity(clu_eigen),4)))
		if((length(table(clu_eigen$membership))>1) & (length(table(clu_eigen$membership))/clu_eigen$vcount!=1)){dendPlot(clu_eigen, mode="hclust")}
	dev.off()
}
# community_detection(net,"full")
community_detection_for_large_graph <- function(in_net, pre){ # This variant of the function community_detection was created for use when a network is too large for the computer's memory, which results in segmentation faults (core dump). This has no edge betweenness calculartion as it is a very expensive algorithm computationally
# 	lo <- layout_with_fr(in_net) # uncomment if a new plotting scheme is required
	write.table(as_long_data_frame(in_net),paste0(prefix,"-cut-",min_val,"-",pre,"-edge_specs.tsv"), sep="\t", quote=FALSE, col.names=NA) # Print the edge table (comment this if required as this way too large for big networks)
	deg <- degree(in_net) # recalculate degree
	clu_louv <- cluster_louvain(in_net) # cluster with louvain algorithm requires only positive values for edge weights
	clu_lab_prop <- cluster_label_prop(in_net) # clustering based on propagating labels (most common edge paths)
	clu_fast_greed <- cluster_fast_greedy(as.undirected(in_net)) # greedy clustering alogorithm
	clu_walktrap <- walktrap.community(in_net,steps=5) # the original article recommends 5 steps
	clu_eigen <- leading.eigenvector.community(in_net)
	# Append edge betweennes
	clusters <- as.data.frame(cbind("louvain"=membership(clu_louv), "lab_prop"=membership(clu_lab_prop), "greedy"=membership(clu_fast_greed), "walktrap"=membership(clu_walktrap), "eigen"=membership(clu_eigen)))
	# Create clusters table
	clusters[nrow(clusters)+1,] <- c(modularity(clu_louv), modularity(clu_lab_prop), modularity(clu_fast_greed), modularity(clu_walktrap), modularity(clu_eigen))
	rownames(clusters)[nrow(clusters)] <- "modularity"
	if(print_names==FALSE){ # If the user requested not to use the original names
		clusters <- as.data.frame(cbind(clusters,"vertex_name"=new_names[rownames(clusters)])) # this uses a simple trick to keep the names of df "clusters"
		clusters <- cbind(clusters[ncol(clusters)],clusters[-ncol(clusters)]) # now shuffle columns into the correct order
	}
	write.table(clusters,paste0(prefix,"-cut-",min_val,"-",pre,"-clusters.tsv"), sep="\t", quote=FALSE, col.names=NA) # 
	pdf(paste0(prefix,"-cut-",min_val,"-",pre,"-clust.pdf"))
		plot(clu_louv, in_net, vertex.size=sqrt(deg),edge.width=E(in_net)$weight, main="Network clusters based on louvain algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)#, edge.width=rescale_vect(E(in_net)$weight,0.5,3)) #, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_louv$membership)), " / Modularity score: ", round(modularity(clu_louv),4)))
		plot(clu_lab_prop, in_net, vertex.size=sqrt(deg), main="Network clustering based on propagating labels algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)#, edge.width=rescale_vect(E(in_net)$weight,0.5,3))#, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_lab_prop$membership)), " / Modularity score: ", round(modularity(clu_lab_prop),4)))
		plot(clu_fast_greed, as.undirected(in_net), vertex.size=sqrt(deg), main="Network clustering based on fast greedy algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo) # , edge.width=rescale_vect(E(in_net)$weight,0.5,3))#, vertex.label=NA)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_fast_greed$membership)), " / Modularity score: ", round(modularity(clu_fast_greed),4)))
		if((length(table(clu_fast_greed$membership))>1) & (length(table(clu_fast_greed$membership))/clu_fast_greed$vcount!=1)){dendPlot(clu_fast_greed, mode="hclust")}
		plot(clu_walktrap, in_net, vertex.size=sqrt(deg), main="Network clustering based on walktrap  algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_walktrap$membership)), " / Modularity score: ", round(modularity(clu_walktrap),4)))
		if((length(table(clu_walktrap$membership))>1) & (length(table(clu_walktrap$membership))/clu_walktrap$vcount!=1)){dendPlot(clu_walktrap, mode="hclust")}
		plot(clu_eigen, in_net, vertex.size=sqrt(deg), main="Network clustering based on leading eigenvector algorithm", vertex.label.color="black", vertex.label.family="Helvetica", layout=lo)
		mtext(paste("Cutoff value: ", min_val, "/ Total Communities: ", length(table(clu_eigen$membership)), " / Modularity score: ", round(modularity(clu_eigen),4)))
		if((length(table(clu_eigen$membership))>1) & (length(table(clu_eigen$membership))/clu_eigen$vcount!=1)){dendPlot(clu_eigen, mode="hclust")}
	dev.off()
}

### MAIN ###
# edge_width(df)
# 1.- Load table and create raw edge map:
edges <- load_edge_map(df) # create a list of all edges
# print(paste("Raw range of values in matrix: ",min(edges[,3]), "to", max(edges[,3]))) # print some useful information for the user to define the cutoff 
# print(paste("Raw Mean:",mean(edges[,3])))
# print(paste("Raw Median:",median(edges[,3])))
# print("Raw quantiles:")
# print(quantile(edges[,3], seq(0,1,0.05)))

# 2.- Sparsification using filters:
edges <- filter_map(edges, is_corr, min_val) # now remove those below the threshold (user-defined). If the user specified this is a correlation matrix, then filter considering absolute values instead. Filtering is considered a sparsification step
# IMPORTANT NOTE: If no sparsification is carried out, most distance-based matrices will create a dense network where most nodes are connected. In these types of networks, weights are more important.

# 3.- Create network
net <- graph_from_data_frame(edges,directed=F, vertices=cbind("names"=rownames(mdf),mdf)) # Now create the network (there are plenty of ways but one of the most controlled ones is to input a table with input node, out node and weight (plus any other metadata for edges. IMPORTANT: names should be added as the first column (which can be easily done with cbind).
# The resulting network should be undirected, weighted (UNW-)
# Now, drop unconnected nodes (conectivity=0) as they affect calculations (disjointed groups may also do this but not as much)
net <- prune_net(net)
lo <- layout_with_fr(net) # Create a global variable for general layout parameters

# 4.- Network analysis
node_meta <- structure_analysis(net, "full") # node_meta creates a list of wo items: the node centrality global table and the global specs, both as df
plot_meta(net, mdf, "full")
# 5.- Community detection
# community_detection(net,"full")
community_detection_for_large_graph(net,"full")

#OTHER tests (DEPRECATED)
# # set.seed(98741987);
# edges2 <- load_edge_map(df2)
# edges2 <- filter_map(edges2, TRUE, 0.2)
# net2 <- graph_from_data_frame(edges2,directed=F, vertices=cbind("names"=rownames(mdf),mdf)) # Now create
# net2 <- prune_net(net2)
# structure_analysis(net2, "full2")
# 
# # , vertex.label=NA
# 
