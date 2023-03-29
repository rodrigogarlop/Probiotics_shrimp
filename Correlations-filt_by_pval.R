# Started 2020-02-11 by Rodrigo García-López
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock" 
# It is intended to create two correlation tables based on a contingency matrix as loades from stdin, once per rows, and once per columns
# The resulting matrix is then filtered by the associated pvalues (cutoff is defined by the user)
# The user must specify whether pearson or spearman should be used
# Please run as follows: cat input_table.tsv|Rscript Correlations.R <out_prefix> <corr_method> <rows_or_cols> <p-value_cutoff>

# Tested with:
# df = read.table("02_split_tables/genetics_lvl2-I_grp_A.tsv",comment.char = "", sep ="\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE)
# out_name = "test"
# method = "pearson"
# rows_or_cols = "rows"
# pval_cut = 0.05

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<4) { # the following arguments are required
  stop("A minimum of 4 arguments are required. Please run as follows: cat input_table.tsv|Rscript Correlations.R <out_prefix> <corr_method> <rows_or_cols> <p-value_cutoff>", call.=FALSE)
}
out_name <- as.character(args[1]) # gets the output prefix (a path may be included)
method <- as.character(args[2]) # gets the correlation method, may be pearson, spearman or kendall
rows_or_cols <- as.character(args[3]) # gets whether rows or cols
pval_cut <- as.numeric(args[4]) # gets whether rows or cols

print("Loaded parameters:")
print(paste("Outfile path/prefix:",out_name))
print(paste("Method:",method))
print(paste("Use:",rows_or_cols))
print(paste("p-value cutoff:",pval_cut))

### Load infile and check parameters ###
df <- read.table(file('stdin'),comment.char = "", sep ="\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE)
opt1 <- c("spearman"=1,"pearson"=1)
opt2 <- c("rows"=1,"cols"=1)
if(is.na(opt1[method])){stop(paste0("Aborting: Please specify the method: pearson, spearman or kendall"), call.=FALSE)}
if(is.na(opt2[rows_or_cols])){stop(paste0("Aborting: Please specify whether columns or rows should be used"), call.=FALSE)}

### Load functions ###
pval_cor <- function(x, y, meth){ # mini-function for getting the pvalue only, using the desire corr method
	pval <- cor.test(x, y, method=meth)$p.value # uses the defined method
	return(pval)
}
cor.test.matrix <- function(mat, meth){ # read a matrix and carry out all permutations to fill a square sym matrix like in cor() function. depends on pval_cor() function above. This uses rows by default
    out <- outer(colnames(mat), colnames(mat), Vectorize(function(i,j) pval_cor(as.numeric(mat[,i]), as.numeric(mat[,j]), meth))) # use outer product (matrix output) to apply the function for pval calculation
    dimnames(out) <- list(colnames(mat), colnames(mat)) # set the names for the square matrix
    return(out)
}

### MAIN ###
if(rows_or_cols=="rows"){
	out_mat <- cor(t(df), method = method) # create the correlation matrix for rows
	pval_mat <- cor.test.matrix(t(df), method)
} else if(rows_or_cols=="cols"){
	out_mat <- cor(df, method = method) # create the correlation matrix for columns
	pval_mat <- cor.test.matrix(df, method)
}
# Optionally, print both the correlations and the pvalues matrices (uncomment if requiered)
# write.table(round(out_mat,2), paste0(out_name,"-corr_",method,"_",rows_or_cols,".tsv"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
# write.table(round(pval_mat,4), paste0(out_name,"-corr-pval_",method,"_",rows_or_cols,".tsv"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

# Finally, use the matching pvalues as a mask for the original output matrix, where those > cutoff should be ingnored (turned to 0)
out_mat <- out_mat*((pval_mat <= pval_cut)*1)

write.table(round(out_mat,2), paste0(out_name,"-corr_",method,"_",rows_or_cols,"_maxpval_",pval_cut,".tsv"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
