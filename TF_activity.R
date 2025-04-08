# TF activity heatmp on control, Full CAR + CAR-E + wt IL2 , cd3only , 4-1bb, icd del
#setwd("/Volumes/My_Passport/paper2_result/new_pipeline/")
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

dir_name<-"TF_CAR_constucts"
samples<- list("Full_CD8_V7", "Costimonly_CD8_V7", "CD3only_CD8_V7", "fulldel_CD8_V7", "Full_CD8_control")

#samples_ind <- lapply(samples, function(sample) grep(sample, colnames(data)))
samples_ind <-lapply(samples, function(sample) grep(sample, colnames(data)))

#all_indices<- unlist(samples_ind)
data_files<- data[, unlist(samples_ind)]

meta<- matrix("", nrow = ncol(data_files) , ncol=1)
colnames(meta)<- c("sample_type" )
meta<- as.data.frame(meta)
sample_names<- colnames(data_files)

rownames(meta)<- sample_names
meta<- as.data.frame(meta)
meta[,1]<- substr(rownames(meta) , 1 , nchar(rownames(meta))-2 )
# meta$samplenames<- c(paste0(group1, "_" , 1:length(Trx)) , paste0(group2, "_" , 1:length(Ctrl))) ## add to code
meta$samplenames<- rownames(meta)
#meta$Batch<- c( rep(batches$Full_CD8_V7 , 4) , rep(batches$dualCAR_CD8_V7 , 4) )
meta

data_files<-data_files%>%rownames_to_column(var="SYMBOL")
rownames(data_files)<- data_files$SYMBOL

data_files$SYMBOL<-NULL
# deseq with rlog 
dds <- DESeqDataSetFromMatrix(data_files, meta, design = ~ sample_type)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

sample_acts <- run_ulm(mat=normalized_counts, net=net, .source='source', .target='target',
                       .mor='mor')
n_tfs <- 40

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
tfs<- c(tfs,"STAT5")
tfs<- unique(tfs)
sample_acts_mat <- sample_acts_mat[,tfs]

sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot

pdf(file.path(  paste0("TF_activity_CD4vsCD8", ".pdf")),height = 10, width = 20)
print(pheatmap(sample_acts_mat, border_color = NA, color=my_color , fontsize_row = 20 , breaks =my_breaks ,fontsize_col = 15 , fontsize = 15))
dev.off()

