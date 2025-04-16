setwd("/Volumes/My_Passport/paper2_result/Github/CAR-E_2/")
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(purrr)
library(DESeq2)
rm (list = ls())

### USER only changes these lines###
dir_name<-"TF_CAR_constructs_V7" # put the name of contrast or figure
samples<- list("Full_CD8_V7","Costimonly_V7","CD3_only_V7","fulldel_CD8_V7") # put conditions of interest
n_tfs <- 50 # number of top most variable TFs among samples
#### END ####



#net <- get_collectri(organism='human', split_complexes=FALSE)
net<- read.table("data/collectRi_net.txt" , header = T)
net$source[net$source=="STAT5A"]<-"STAT5"
net$source[net$source=="STAT5B"]<-"STAT5"
net$target[net$target=="STAT5A"]<-"STAT5"
net$target[net$target=="STAT5B"]<-"STAT5"
net$source[net$source=="REL"]<-"RELA"
net$target[net$target=="REL"]<-"RELA"

sub_net<- net[,1:2]
ind<- which(duplicated(sub_net))
net<- net[-ind, ]

source("create_conditions.R")
data<- read.table("data/countMatrix_all.txt" , header = T)
geneIDs<- data$Geneid
gene_IDs <- read.table("data/Ensemble_IDs.txt", header = T , fill = T)
gene_IDs[gene_IDs==""]<-NA
colnames(gene_IDs) <- c("ensembl_gene_id", "hgnc_symbol")

data$Geneid<- NULL

rownames(data)<- geneIDs
data<- data[complete.cases(data),]
geneIDs<- rownames(data)

data <- data %>%
  rownames_to_column(var = "ensembl_gene_id")
data <- data %>%
  left_join(gene_IDs, by = "ensembl_gene_id")
data$hgnc_symbol[data$hgnc_symbol==""]<-NA
data<- data[complete.cases(data), ]

ind.dup<- which(duplicated(data$hgnc_symbol))
if (length(ind.dup)>0)
{
  data<- data[-ind.dup,]
}

Genes<- data$hgnc_symbol
ENSEMBLE<-data$ensembl_gene_id
data$ensembl_gene_id<-NULL
rownames(data)<- Genes
data$hgnc_symbol<-NULL

# get normalized counts

conditions<- creat_conditions(data)


if (!dir.exists("TF_activity"))
  dir.create("TF_acti")


samples_ind <-lapply(samples, function(sample) grep(sample, colnames(data)))
data_files<- data[, unlist(samples_ind)]

meta<- matrix("", nrow = ncol(data_files) , ncol=1)
colnames(meta)<- c("sample_type" )
meta<- as.data.frame(meta)
sample_names<- colnames(data_files)

rownames(meta)<- sample_names
meta<- as.data.frame(meta)
meta[,1]<- substr(rownames(meta) , 1 , nchar(rownames(meta))-2 )
meta$samplenames<- rownames(meta)
meta

# data_files<-data_files%>%rownames_to_column(var="SYMBOL")
# rownames(data_files)<- data_files$SYMBOL

#data_files$SYMBOL<-NULL
# deseq with rlog 
dds <- DESeqDataSetFromMatrix(data_files, meta, design = ~ sample_type)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

sample_acts <- run_ulm(mat=normalized_counts, net=net, .source='source', .target='target',
                       .mor='mor')


# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# V7_ctrl_select<- c(paste0("Full_CD8_V7_", 1:4), paste0("Full_CD8_control_", 1:4))
# tfs <- sample_acts[sample_acts$condition %in% V7_ctrl_select,] %>%
#   group_by(source) %>%
#   summarise(std = sd(score)) %>%
#   arrange(-abs(std)) %>%
#   head(n_tfs) %>%
#   pull(source)


tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
tfs<- c(tfs,"STAT5")
tfs<- unique(tfs)
ind<- which(colnames(sample_acts_mat)%in% tfs)
sample_acts_mat <- sample_acts_mat[,ind]

sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

anno_col<- sample_acts_mat%>%t()%>% colnames()%>% substr(. , 1, nchar(.)-2)%>% as.data.frame()
rownames(anno_col)<- rownames(sample_acts_mat)
# Plot

pdf(file.path(  paste0("TF_activity/",dir_name, ".pdf")),height = 20, width = 10)

print(pheatmap(t(sample_acts_mat), border_color = NA, color=my_color , fontsize_row = 20 ,annotation_col =  anno_col,
               breaks =my_breaks ,fontsize_col = 20 , fontsize = 15, labels_col = NULL, labels_row = NULL))
dev.off()

