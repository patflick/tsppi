# TODO: Add comment
# 
# Author: flick
###############################################################################

if (.Platform$OS.type == "unix")
{
	setwd("~/dev/ppi/hpa");
} else {
	setwd("Y:\\dev\\ppi\\hpa");
}

#file = "data/normal-tissue.csv"
cerebral = "cerebral.cortex_neuronal.cells.csv"
liver = "liver_hepatocytes.csv"


fpath = file.path(getwd(), "data", "tissues", cerebral);
normal_tissue_data <- read.csv(fpath, header=TRUE, sep=",");

cerebral_which <- which(normal_tissue_data$Expression.type == "APE" & (normal_tissue_data$Reliability == "High" | normal_tissue_data$Reliability == "Medium"))
cerebral_genes <- normal_tissue_data$Gene[cerebral_which]


fpath = file.path(getwd(), "data", "tissues", liver);
liver_data <- read.csv(fpath, header=TRUE, sep=",");

liver_which <- which(liver_data$Expression.type == "APE" & (liver_data$Reliability == "High" | liver_data$Reliability == "Medium"))
liver_genes <- liver_data$Gene[liver_which]

