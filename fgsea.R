#Reference: https://biostatsquid.com/fgsea-tutorial-gsea/

set.seed(349)
library(fgsea)
# GSEA for many genesets  -------------------------
in_path <- "mSigDB/"
out_path <- "fGSEA results/"
# Functions Background Gene Prep ===================================================
## Function: Adjacency matrix to list 
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt 
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('.gmt conversion successfull')
  return(final_list)
}



## 1. Read in data 
list.files(in_path)
df <- read.csv()
names(df)[names(df) == "V2"] <- "gene_symbol"
names(df)[names(df) == "V1"] <- "Ensembl ID"
## 2. Prepare background genes 
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene_symbol
list.files(in_path)
gmt_files <- list.files(path = in_path, pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files[1], my_genes, savefile = FALSE) #PMID: 26346307
## 3. Rank Genes
# RANK GENES 
rankings <- sign(df$avg_log2FC)*(-log10(df$p_val)) 
names(rankings) <- df$gene_symbol 
rankings <- sort(rankings, decreasing = T)
plot(rankings)
## Fix Infinite values 
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)
# Save ranking file
write.table()

## 1. Read in data 
list.files(in_path)
df <- read.csv("with ensemble/Fibroblast_tumors_KPCSvKPC_2025.05.23_ensembl.csv")
names(df)[names(df) == "V2"] <- "gene_symbol"
names(df)[names(df) == "V1"] <- "Ensembl ID"
## 2. Prepare background genes 
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene_symbol
list.files(in_path)
gmt_files <- list.files(path = in_path, pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files[1], my_genes, savefile = FALSE) #PMID: 26346307
## 3. Rank Genes
# RANK GENES 
rankings <- sign(df$avg_log2FC)*(-log10(df$p_val)) 
names(rankings) <- df$gene_symbol 
rankings <- sort(rankings, decreasing = T)
plot(rankings)
## Fix Infinite values 
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)
# Save ranking file as .rnk
write.table()
  
##4. Run GSEA
x <- #.gmt file of interest
bg_genes <- prepare_gmt(gmt_files[1], my_genes, savefile = TRUE)
writeGmtPathways(bg_genes, x)
GSEA <- fgsea(pathways = bg_genes,
                stats = rankings,
                scoreType = 'std',
                minSize = 5,
                maxSize = 500)
  ##name of comparison
  name_of_comparison <- ''
  background_genes <- ''
  filename <- paste0(out_path, name_of_comparison, '_', background_genes)
  fwrite(GSEA, file = paste0(filename, '_gsea_results.tsv'))
  GSEAcsv <- apply(GSEA,2,as.character)
  write.csv(GSEAcsv, file = paste0(filename, '_gsea_results.csv'))

