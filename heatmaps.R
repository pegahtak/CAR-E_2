rm(list=ls())
#setwd("/Volumes/My_Passport/paper2_result/new_pipeline/DEG/")
library(dplyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(purrr)
library(readxl)
library(tibble)
library(xlsx)
library(RColorBrewer)
library(tibble)
rm(list = ls())

# define contrasts samples and control
main<- "Full_CD8_V7"
control<- "fulldel_CD8_V7"
other_cond<- c("Costimonly_CD8_V7", "CD3only_CD8_V7")
contrast<- paste0(main, "_vs_", control)
n_top<- 40
fig_name<- "6-c"
PATH="/Volumes/My_Passport/paper2_result/new_pipeline/DEG/" # put current path

# read DEGs in main vs control 
geneDat<- read.xlsx(paste0(contrast,"/DEGs_", contrast, ".xlsx"), sheetIndex=1)
geneDat$NA.<-NULL

#select top genes
DEGs<- geneDat %>% arrange(desc(log2FoldChange) )%>%pull(Symbol)
topGenes <- geneDat %>% arrange(desc(log2FoldChange) )%>%pull(Symbol)%>% head(n_top+10)


sample_names<- rev(c(main, other_cond))# using rev to put the first column on right


files_names<- paste0(sample_names, "_vs_",control )
list_files<-as.list(paste0(PATH , files_names, "/res_lfc_", files_names , ".txt"))
x_list <- lapply(list_files, read.delim)

# Make a data that contians LFC of the selected top genes in all conditions vs control
x <- lapply(x_list, function(file) {
  subset(file, Symbol %in% topGenes)[, c("log2FoldChange", "Symbol"), drop = FALSE]
})
x<- lapply(x, function(file) file[!duplicated(file$Symbol), ])

final_table <- do.call(cbind, x)
typeof(final_table)
rownames(final_table)<- NULL
final_table<- column_to_rownames(final_table,var="Symbol")
final_table<- final_table[, grep("log2FoldChange", colnames(final_table)) ]


final_table<- final_table %>%
  mutate(sd = rowSds(as.matrix(across(everything()))))%>%
  arrange(desc(sd))%>% as.data.frame()%>%dplyr::select(-sd)

#topGenes250<- rownames(final_table)%>%head(250)
#topGenes50<- final_table%>%rownames()%>%head(50)

#rownames(final_table)<- final_table$Symbol
#Symbols<-  final_table$Symbol
#final_table<- final_table[,-grep("Symbol" , colnames(final_table))]
colnames(final_table)<- files_names

# remove rows with NA
#final_table <- final_table[complete.cases(final_table), ]
#write.csv(final_table, paste0(fig_name,".csv")  , quote = F)

#final_table<- final_table[ order(final_table[,ncol(final_table)] , decreasing = T), ]
colnames(final_table)<- sample_names
final_table<- final_table[which(rownames(final_table)%in%topGenes),]
final_table<- final_table[complete.cases(final_table),]

#final_table<- final_table%>% arrange(desc(!!sym(main) ))%>% head(n_top)
t<- final_table
if (! dir.exists("heatmaps"))
  dir.create("heatmaps")

write.csv(final_table ,paste0("heatmaps/", fig_name, ".csv") , quote = F)


# removing outliers
t[t > 20] <- NA
#t[t< -12]<- 
max<- max(t)
t[t==NA]<- max
t<- t[complete.cases(t),]
t<- t%>% arrange(desc(!!sym(main) ))%>% head(n_top)


# r_ind<- which(rownames(final_table)%in%remove)
# t<- final_table[-r_ind, ]

mycols<- colorRampPalette(c("blue", "white", "red"))(100)
# myBreaks<- c(-4, -2, 4, 10)
#ph.breaks <- seq(from = -10, to = 10, length.out = 201)
t_mean<- t%>% unlist()%>% mean()
ph.breaks <- c(seq(from = min(t, na.rm = TRUE), to = t_mean , length.out = 50), 
               seq(from = t_mean+0.01, to = max(t, na.rm = TRUE) + 0.01, length.out = 51))
pdf(paste0("heatmaps/", fig_name , ".pdf"), width = 4 , height = 12)
#par(mar=c(20,20,20,20))
print(pheatmap(t , cluster_cols = F , show_colnames = T,
               fontsize_row = 15 , fontsize_col = 15 , cluster_rows = T, fontsize = 15 , legend = T
               , annotation_legend = T , angle_col = 315 , na_col="purple" , clustering_method = "average", breaks = ph.breaks ))
dev.off()
graphics.off()