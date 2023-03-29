# Started on 2020-09-22
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# The script is intended to extract groups a contingency table based on the sample names
# Groups are defined by the user
# The input is contingency table with feature names in column 1 and sample names in row 1

# Run as follows:
# cat table.tsv|Rscript Extract_groups_from_table <prefix_output> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>

# Tested with command:
# cat Filter_core_H/rarefied_7k-3y-H-lvl6.tsv|Rscript Rscript Extract_groups_from_table test HE IE HL IL HM IM
# Test in R:
# df <- read.table("Filter_core_H/rarefied_7k-3y-H-lvl6.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
# prefix="test"
# groups <- c("Z1","Z2","Z3","Z4","Z5")
# groups <- c("2015","2016","2017","2018")
# groups <- c("H","I","S")
# groups <- c("H")
# groups <- c("HE","IE","HL","IL","HM","IM")
# groups <- c("E","L","M")
# groups <- c("2015_Z1_H","2015_Z1_I","2015_Z2_H","2015_Z2_I","2016_Z2_H","2016_Z2_I","2015_Z3_I","2015_Z3_S","2016_Z3_S","2017_Z3_H","2017_Z3_I","2017_Z3_S","2018_Z4_I","2018_Z5_H")
# # groups <-("2015_Z1_I", "2015_Z2_H", "2016_Z2_H", "2016_Z2_I", "2017_Z3_H", "2017_Z3_I", "2018_Z4_I", "2018_Z5_H")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) { # at least, two arguments should be included: <min_nonNAs> <prefix_output>  <name_of_alpha_metric>
  stop("A minimum of 2 arguments are mandatory: cat table.tsv|Rscript PCoA_from_dm.R <prefix_output>  <name_of_beta_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
groups <- as.character(args[2:length(args)]) # Create a vector with the name of the groups (this should be included in the sample names and be exclusive to each group for this to work)
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # 
print("Loaded parameters:")
print(paste("Prefix:",prefix))
print("Groups that should be present:")
print(groups)

####################################### Pre-load #######################################
### Define Functions ###
# Predefine some important global objects
prepare_data <- function(df) { # create presets (sample totals and groups # MOD: ignored empty groups
	samples <- names(df) # extract names
	grp_sam <- sapply(groups,function(x) grep(x,samples)) # create index of samples in each group
	grp_less <- sapply(grp_sam,length) # and count the totals # This will now consider empty groups
	grp_all <- grp_less # Since this was modified, we still required the old version
	grp_sam <- grp_sam[(grp_less>=1)] # Remove groups that ended up with 0 samples
	grp_less <- sapply(grp_sam,length) # This will now consider only non-empty groups
	len <- length(unlist(grp_sam)) # Get total items
	grp <- rep("Other",len) # create a template for considering non-group items
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual sample distribution for all groups
		for(j in 1:length(grp_sam[[i]])){
		grp[grp_sam[[i]][j]] <- names(grp_less)[i];
		}
	}
	out <-list(grp_sam,grp_less,as.factor(grp),grp_all)
	names(out) <- c("grp_sam","grp_less","grp", "grp_all")
	return(out)
}

####################################### MAIN #######################################
############ Part 1: Assess the group distribution, calculate stats ############
 ### Prefilter samples ###
if(length(groups)<2){groups <- c(groups,"not_found")} # Patch with whatever is not found in the sample names to fix the one group case
presets <- prepare_data(df) # Define how samples from each group are distributed
# Remove samples not in groups requested by the user:
df <- df[,unlist(presets$grp_sam)] # filter columnwise
df <- df[rowSums(df)>0,]# filter columnwise
write.table(df,paste0(prefix,"_grp_",gsub(", ","_",toString(names(presets$grp_less))),".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) 
