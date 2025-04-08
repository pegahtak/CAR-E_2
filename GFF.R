library(rtracklayer)
system("wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz -P data/")
system("gunzip data/Homo_sapiens.GRCh38.113.chr.gtf.gz ")
gff<- readGFF("data/Homo_sapiens.GRCh38.113.chr.gtf")
p_genes<- gff[gff$gene_biotype=="protein_coding",c("gene_id", "gene_name") ]
system("rm data/Homo_sapiens.GRCh38.113.chr.gtf")
# remove duplicates
p_genes<- p_genes[!duplicated(p_genes), ]
rm (gff)
