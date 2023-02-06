# Capture Command Line Arguments
cmdArgs <- commandArgs(trailingOnly = TRUE)
filesfilter <- cmdArgs[2]
filesgene <- cmdArgs[3]
filesoutput <- cmdArgs[4]

## filter 400kb window
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(coloc)

setwd(filesfilter)
my_files <- fread(file = filesgene, header = FALSE)
my_data <- list()
for (i in seq_along(my_files$V1)) {
	my_data[[i]] <- fread(file = my_files$V1[i], header =FALSE)
	my_data[[i]] <- subset(my_data[[i]], my_data[[i]]$V3 < 200000 & my_data[[i]]$V3 > -200000)
}
setwd(filesoutput)
for (i in seq_along(my_files$V1)) {
	lapply(my_data[[i]], function(x){write.table(my_data[[i]]$V1, file = paste(my_files[i], ".400kbtxt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep= "\t")})
}
