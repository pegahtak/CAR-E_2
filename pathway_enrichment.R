setwd("/Volumes/My_Passport/paper2_result/Github/CAR-E_2/")
library (org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(tidyverse)
library(DOSE)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(pathview)
library(genekitr)
library(signatureSearch)
library(enrichplot)


rm(list = ls())
#setwd("/Volumes/My_Passport/paper2_result/new_pipeline/DEG/")
source("GFF.R") # creates p_genes with ensemble IDs and gene symbols
# read gene sets
c.hall <- read.gmt("data/h.all.v2024.1.Hs.symbols.gmt")
c.7<- read.gmt("data/c7.all.v2024.1.Hs.symbols.gmt")

source("create_conditions.R")
data<- read.table("data/countMatrix_all.txt", header = T)
geneIDs<- data$Geneid
gene_IDs <- read.table("/Volumes/My_Passport/paper2_result/new_pipeline/Ensemble_IDs.txt", header = T , fill = T)
gene_IDs[gene_IDs==""]<-NA
colnames(gene_IDs) <- c("ensembl_gene_id", "hgnc_symbol")

data$Geneid<- NULL

conditions<-creat_conditions(data)

rownames(data)<- geneIDs
colnames(data)<- unlist(map(strsplit(colnames(data), split = "_CKD"), 1))
data<- data[complete.cases(data),]
geneIDs<- rownames(data)



#setwd("/Volumes/My_Passport/paper2_result/new_pipeline/pathway_analysis/")
getwd()

for ( i in 1:(length(conditions) - 1))
{
  for( j in (i+1):(length(conditions)))
  {
    # get indices
    if ( i == j)
      next
    
    Trx<- conditions[[i]]
    Ctrl<- conditions[[j]]
    
    group1<- names(conditions)[i]
    group2<- names(conditions)[j]
    
    dir_name<- paste0(group1 , "_vs_" , group2)
    #dir.create(paste0("pathway_analysis/",dir_name) , showWarnings = T)
    
    # read deseq2 results 
    res<- read.delim(paste0( group1, "_vs_", group2,"/res_lfc_",group1, "_vs_", group2 , ".txt" )
                     , header = T)
    
    # remove non coding genes or genes with NA values in padj.
    res<- res[res$Symbol%in%p_genes$gene_name & ! is.na(res$padj),]
    
    
    ids <- bitr(res$Symbol, 
                fromType = "SYMBOL", 
                toType = c("ENSEMBL", "ENTREZID"), 
                OrgDb = "org.Hs.eg.db")
    
    
    # add entrez IDs and remove rows with missing pvalue
    x<-res%>%filter(Symbol%in%ids$SYMBOL)%>%left_join(ids, by = c("Symbol" = "SYMBOL"))%>%
      filter(!is.na(pvalue))
    
    # remove duplicates
    x<-x[which(duplicated(x$ENTREZID) == F), ]
    foldchanges<- x$log2FoldChange
    names(foldchanges)<- x$ENTREZID
    foldchanges<- sort(foldchanges , decreasing = T)
    
    msig.hall <- GSEA(foldchanges, TERM2GENE=c.hall, verbose=FALSE, minGSSize=15,maxGSSize=500
                      ,pvalueCutoff=0.05 )
    msig_df_hall <- data.frame(msig.hall)
    msig_df_hall<- msig_df_hall[order(msig_df_hall$NES , decreasing = T), ]
    write.xlsx(msig_df_hall,  file.path(paste0(dir_name ,paste0("/",dir_name, "_C7_MsigDB.xlsx")) ) )
    
    pdf(paste0(dir_name,"/IL2_STAT5_SIGNALING.pdf") )
    print(gseaplot(msig.hall , geneSetID = "HALLMARK_IL2_STAT5_SIGNALING"))
    dev.off()
    
  }
}

