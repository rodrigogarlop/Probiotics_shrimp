# This will be a fixed script that will only work for plotting the density plot of the distributions comparing A and L ponds for each organ.
# We want to compare all centrality measurements at each tax lvl
# It won't be usefull for any other comparison.
library("ggplot2")
library("ggpubr")
setwd("/home/rod/Documents/02_Collaborations/Prebioticos_Pablo/cor_0.8_pval_0.01/")
lvl <- c("OTU")
# tab1 <- "04_networks/prebiotics_lvl5-H_grp_A-corr_pearson_rows-cut-0.5-full-node_specs.tsv"
# tab2 <- "04_networks/prebiotics_lvl5-H_grp_L-corr_pearson_rows-cut-0.5-full-node_specs.tsv"
compare_tables <- function(tab1, tab2, tab3, column){ # Use this to go through each three files
	df1 <- read.table(tab1, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # First, load each file
	df2 <- read.table(tab2, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
	df3 <- read.table(tab3, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
	df1 <- df1[1:(nrow(df1)-2),] # and remove global features
	df2 <- df2[1:(nrow(df2)-2),]
	df3 <- df3[1:(nrow(df3)-2),]
	rownames(df1) <- df1[,1] # also, vertices names (v1...n don't mean anything), so replace them with actual taxa/cluster names
	rownames(df2) <- df2[,1]
	rownames(df3) <- df3[,1]
	temp <- table(c(rownames(df1),rownames(df2),rownames(df3))) # Now create a single vector
	temp <- names(temp[temp==3]) # to determine which elements are present in all three tables
	df1 <- df1[temp,] # Now, filter matching items in all three sets
	df2 <- df2[temp,]
	df3 <- df3[temp,]
	# Also, remove rows where NAs are present (missing value in any group) as they cannot be compared
	out <- cbind(df1[column], df2[column], df3[column]) # Extract the target column as data.frame
	out <- out[!is.na(rowSums(out)),] # remove rows having at least one NA
	return(out)
}
trans_col <- function(color, percent = 50, name = NULL) { # Set transparent color vector: Function from ## www.dataanalytics.org.uk
	rgb.val <- col2rgb(color)
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100 - percent) * 255 / 100, names = name)
	invisible(t.col) ## Save the color
}
plot_density <- function(tab, lvl){ # from a 2-column table, plot densities
	a <- tab[,1][!is.na(tab[,1])]
	b <- tab[,2][!is.na(tab[,2])]
	c <- tab[,3][!is.na(tab[,3])]
	ac <- ((!is.na(tab[,1])))*((!is.na(tab[,3]))) # determine which are directly comparable (ocurring in 2 gprs)
	bc <- ((!is.na(tab[,2])))*((!is.na(tab[,3])))
	A <- stats::density(a)
	B <- stats::density(b)
	C <- stats::density(c)
	range.x <- range(c(A[[1]],B[[1]]),C[[1]])
	max.y <- max(c(A[[2]],B[[2]]),C[[2]])
	wpval.ac <- wilcox.test(tab[as.logical(ac),1],tab[as.logical(ac),3], paired=T)$p.value
	tpval.ac <- t.test(tab[as.logical(ac),1],tab[as.logical(ac),3], paired=T)$p.value
	cor.ac <- round(cor(tab[as.logical(ac),1],tab[as.logical(ac),3]),4)
	if(wpval.ac < 0.0001){wpval.ac = "< 0.0001"}else{wpval.ac=round(wpval.ac,4)}
	if(tpval.ac < 0.0001){tpval.ac = "< 0.0001"}else{tpval.ac=round(tpval.ac,4)}
	wpval.bc <- wilcox.test(tab[as.logical(bc),2],tab[as.logical(bc),3], paired=T)$p.value
	tpval.bc <- t.test(tab[as.logical(bc),2],tab[as.logical(bc),3], paired=T)$p.value
	cor.bc <- round(cor(tab[as.logical(bc),2],tab[as.logical(bc),3]),4)
	if(wpval.bc < 0.0001){wpval.bc = "< 0.0001"}else{wpval.bc=round(wpval.bc,4)}
	if(tpval.bc < 0.0001){tpval.bc = "< 0.0001"}else{tpval.bc=round(tpval.bc,4)}
# 	pdf("test.pdf")
	plot(A, xlim=range.x, ylim=c(0,max.y), yaxt='n', xlab="", ylab=paste0("Lvl ",lvl), main="") # Create the density graph
	polygon(A,col=transp_cols[1],border=solid_cols[1]) # and color it accordingly
	polygon(B,col=transp_cols[2],border=solid_cols[2])
	polygon(C,col=transp_cols[3],border=solid_cols[3])
	legend("topleft", legend=(c(paste0("Wilcoxon p-value: ",wpval.ac), paste0("t p-value: ",tpval.ac), paste0("Pearson corr: ",cor.ac))), bty = "n", cex=0.8, title="2% vs C")
	legend("topright", legend=(c(paste0("Wilcoxon p-value: ",wpval.bc), paste0("t p-value: ",tpval.bc), paste0("Pearson corr: ",cor.bc))), bty = "n", cex=0.8, title="10% vs C")
# 	dev.off()
}
plot_ggboxplot <- function(tab) {
	x <- combn(c("Agavin 2%", "Agavin 10%", "Control"),2)
	p_value <- lapply(seq_len(ncol(x)), function(i) x[,i])
	long_format <- stack(tab)
	borders <- c("coral1","cornflowerblue","turquoise1")
	colors <- c("bisque","lightblue","darkslategray1")
	levels(long_format[,2]) <- c("Agavin 2%", "Agavin 10%", "Control")
	boxplot <- ggplot(long_format, aes(x = ind, y = values, color = ind)) + 
		geom_boxplot(fill=colors, color=borders) + theme_classic() +
		theme(legend.position = "none", text = element_text(size = 16), axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 18), plot.title = element_text(hjust = 0.5, size=22)) + 
		labs(title= names(tab)[1], x = "", y = "")  +
		stat_compare_means(comparisons = p_value, method = "wilcox.test") +  stat_compare_means(label.y = max(tab)+(max(tab)-min(tab))*0.3)
	print(boxplot)
}
process_column <- function(column,organ){ # this is just to circle through columns and plot
	for(i in lvl){
		tab1 <- paste0("04_networks/prebiotics_lvl",i,"-",organ,"_grp_2-corr_pearson_rows-cut-0.7-full-node_specs.tsv")
		tab2 <- paste0("04_networks/prebiotics_lvl",i,"-",organ,"_grp_10-corr_pearson_rows-cut-0.7-full-node_specs.tsv")
		tab3 <- paste0("04_networks/prebiotics_lvl",i,"-",organ,"_grp_C-corr_pearson_rows-cut-0.7-full-node_specs.tsv")
		tab <- compare_tables(tab1, tab2, tab3, column)
# 		plot_density(tab,i)
		plot_ggboxplot(tab)
	}
}

 ### Main ###
transp_cols <- c(trans_col("pink",50),trans_col("firebrick",50), trans_col("cornflowerblue",50))
solid_cols <- c("indianred2","red","blue")
pdf("05_collated_comparisons/degree_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(2,"H"); dev.off()
pdf("05_collated_comparisons/degree_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(2,"i"); dev.off()
pdf("05_collated_comparisons/eigenvector_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(3,"H"); dev.off()
pdf("05_collated_comparisons/eigenvector_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(3,"i"); dev.off()
pdf("05_collated_comparisons/closeness_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(4,"H"); dev.off()
pdf("05_collated_comparisons/closeness_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(4,"i"); dev.off()
pdf("05_collated_comparisons/betweenness_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(5,"H"); dev.off()
pdf("05_collated_comparisons/betweenness_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(5,"i"); dev.off()
pdf("05_collated_comparisons/transitivity_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(6,"H"); dev.off()
pdf("05_collated_comparisons/transitivity_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(6,"i"); dev.off()
pdf("05_collated_comparisons/closeness_w-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(7,"H"); dev.off()
pdf("05_collated_comparisons/closeness_w-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(7,"i"); dev.off()
pdf("05_collated_comparisons/betweenness_w-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(8,"H"); dev.off()
pdf("05_collated_comparisons/betweenness_w-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(8,"i"); dev.off()
pdf("05_collated_comparisons/hub_score-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(9,"H"); dev.off()
pdf("05_collated_comparisons/hub_score-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(9,"i"); dev.off()
pdf("05_collated_comparisons/vertex_strength-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(10,"H"); dev.off()
pdf("05_collated_comparisons/vertex_strength-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(10,"i"); dev.off()


# [1] "vertex_name"     "degree_u"        "eigenvector_u"   "closeness_u"    
#  [5] "betweenness_u"   "transitivity_u"  "closeness_w"     "betweenness_w"  
#  [9] "hub.score"       "vertex.strength"

# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
