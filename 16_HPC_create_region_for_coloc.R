library(data.table)
library(rtracklayer)

print(commandArgs(trailingOnly=TRUE))

#Load in the args from command line
args <- commandArgs(trailingOnly=TRUE)

for (arg in args) {
  # skip blank or whitespace-only arguments
  if (nchar(trimws(arg)) == 0) next
  
  ta <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(ta) != 2 || nchar(ta[2]) == 0) {
    stop("Not all arguments are given or properly formatted (key=value)")
  }
  assign(ta[1], ta[2])
}

phenoname <- as.character(phenoname)
input_data_path <- as.character(input_data_path)
output_data_rootname <- as.character(output_data_rootname)
plot_title <- as.character(plot_title)


print("input_data_path:")
print(input_data_path)

full_input_dir <- paste(input_data_path, "/", phenoname, ".assoc.linear", sep = "")

print("full_input_dir:")
print(full_input_dir)

#read in the association data to create the Manhattan plot 
data<-fread(full_input_dir) #plink output format
data<-as.data.frame(data)



top <- data[order(data$P, decreasing = F)[1],]


filt <- data[intersect(which(data$CHR == top$CHR),
                              which(data$BP < (top$BP + 500000) &
                                    data$BP > (top$BP - 500000))),]

filt$BP <- as.numeric(filt$BP)
filt_ordered <- unique(filt[order(filt$BP, decreasing = F),])

write.table(filt_ordered, row.names = F, col.names = T, file= paste0(output_data_rootname, "/",  phenoname, "_for_coloc.txt"))

#note if you need to liftover there is code in script 4