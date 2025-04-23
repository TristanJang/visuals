library(tidyverse)
library(data.table)

read_preprocess_file <- function(isolate) {
  variant_calling_file<-fread(paste("/home/tjang/mnt/lpcdrp/data/depot/assembly/south_africa/", isolate, "/", isolate, ".vcf", sep = ""))
  variant_calling_file$isolate <- isolate
  return(variant_calling_file)
}
isolate <- fread("/home/tjang/mnt/lpcdrp/data/depot/assembly/south_africa/isolates.txt", header = FALSE)
lineage <- fread("/home/tjang/mnt/lpcdrp/data/depot/assembly/south_africa/lineages.tsv")
isolate<- isolate[-c(14,22,3),]
variant_calling_files <- mapply(read_preprocess_file, isolate$V1 , SIMPLIFY = FALSE)
variant_calling_files_unzipped <- do.call(rbind, variant_calling_files)
variant_calling_files_unzipped<- variant_calling_files_unzipped%>%mutate(var_type = if_else(nchar(REF) == nchar(ALT),"snp","indel"))
variant_calling_files_unzipped %>% group_by(isolate,var_type) %>% summarise(number_of_variants=n())%>%
 ggplot() + geom_bar(aes(reorder(isolate,number_of_variants),number_of_variants, fill=var_type),stat = "identity") +
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("isolate name")+ ylab("Number of variants")


read_preprocess_cov_file <- function(isolate) {
  cov_file<-fread(paste("/home/tjang/mnt/lpcdrp/data/depot/assembly/south_africa/", isolate ,
                        "/coverage_depth/",isolate,"_coverage.csv",sep = ""))
  cov_file$isolate <- isolate
  return(cov_file)
}
cov_files <- mapply(read_preprocess_cov_file, isolate$V1, SIMPLIFY = FALSE)
cov_files_unzipped <- do.call(rbind, cov_files)
# generate coverage depth box plots
cov_files_unzipped%>% ggplot(aes(x=V3, color=isolate))+geom_boxplot()+
  xlab("coverage_depth") + ylab("isolate")
# generate coverage summary stats
cov_files_unzipped %>% group_by(isolate) %>% summarise(mean_cov = mean(V3), median_cov = median(V3), sd(V3))


