# Started 2020-02-05 by Rodrigo García-López
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock" 
# It is intended to create two correlation tables based on a contingency matrix as loades from stdin, once per rows, and once per columns
# The user must specify whether pearson or spearman should be used
# Please run as follows: cat input_table.tsv|Rscript Correlations.R <out_prefix> <corr_method> <rows_or_cols>

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<3) { # five arguments are required
  stop("A minimum of 3 arguments are required. Please run as follows: cat input_table.tsv|Rscript Correlations.R <out_prefix> <corr_method> <rows_or_cols>", call.=FALSE)
}
out_name <- as.character(args[1]) # gets the output prefix (a path may be included)
method <- as.character(args[2]) # gets the correlation method, may be pearson, spearman or kendall
rows_or_cols <- as.character(args[3]) # gets whether rows or cols

print("Loaded parameters:")
print(paste("Outfile path/prefix:",out_name))
print(paste("Method:",method))
print(paste("Use:",rows_or_cols))

df <- read.table(file('stdin'),comment.char = "", sep ="\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE)
opt1 <- c("spearman"=1,"pearson"=1)
opt2 <- c("rows"=1,"cols"=1)
if(is.na(opt1[method])){stop(paste0("Aborting: Please specify the method: pearson, spearman or kendall"), call.=FALSE)}
if(is.na(opt2[rows_or_cols])){stop(paste0("Aborting: Please specify whether columns or rows should be used"), call.=FALSE)}
if(rows_or_cols=="rows"){
	out_mat <- cor(t(df), method = method) # create the correlation matrix for rows
}
if(rows_or_cols=="cols"){
	out_mat <- cor(df, method = method) # create the correlation matrix for columns
}
write.table(round(out_mat,2), paste0(out_name,"-corr_",method,"_",rows_or_cols,".tsv"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
