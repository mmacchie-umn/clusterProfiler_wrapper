#!/usr/bin/env Rscript



#######################################################################
# Install and load required R packages
#######################################################################

req_pkgs <- list("tidyverse","clusterProfiler","RColorBrewer","glue","msigdbr","pathview","readr","stringr","enrichplot","ggplot2","ReactomePA")

for (i in req_pkgs){
  if (i %in% rownames(installed.packages())){
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }else{
    BiocManager::install(i,character.only = TRUE)
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }
}



#######################################################################
# Parse options
#######################################################################

# Clear out any objects from workspace
rm(list = ls(all.names = TRUE))


library("optparse",verbose=FALSE,quietly=TRUE)
option_list = list(
  make_option(c("-f", "--file"),
              type="character",
              default=NULL,
              help="Input DE expression result file name (Required)",
              metavar="character"),
  make_option(c("-s","--species"),
              type="character",
              default=NULL,
              help="species name (Required) - underscores should replace spaces. eg. Canis_lupus_familiaris ",
              metavar="character"),
  make_option(c("-d","--dir"),
              type="character",
              default="./",
              help="project directory",metavar="character"),
  make_option(c("-i","--geneIDtype"),
              type="character",
              default="ENSEMBL",
              help="Either ENSEMBL,SYMBOL,ENTREZ...",metavar="character"),
  make_option(c("-n","--num_top_terms"),
              type="numeric",
              default="Inf",
              help="Report top N ontology terms",metavar="character"),
   make_option(c("-g","--numDEGs"),
              type="numeric",
              default="200",
              help="Use top X DE genes for Ontology Term enrichment Analysis. Genes with p-value less than 0.05 (or user-specified cutoff) are sorted (largest to smallest) by the absolute value of log2FC. Top 200 genes are kept by default.",metavar="character"),
  make_option(c("-p","--FDRGenes"),
              type="numeric",
              default="0.05",
              help="FDR cutoff used for differentially expressed genes",metavar="character"),
  make_option(c("-a","--pvalueOnt"),
              type="numeric",
              default="0.05",
              help="p-value cutoff used for ontology term enrichment test",metavar="character"),
  make_option(c("-c","--foldchange"),
              type="numeric",
              default="1",
              help="Log2(fold change) cutoff for differentially expressed genes",metavar="character"),
 make_option(c("-b","--gene_ontology"),
              type="character",
              default="BP",
              help="Biological process (BP), Cellular Component (CC), Molecular function (MF), or all ontologies (ALL)",metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("DE expression results table file is not specified (-f).", call.=FALSE)
}

if (is.null(opt$species)){
  print_help(opt_parser)
  stop("No species is specified (-s). Provide genus species name with underscores. eg. Mus_musculus or Canus_lupus_familiaris.", call.=FALSE)
}


# DEG file 
DEfile <-opt$file
# Number of Ontology terms to report
num_terms <- opt$num_top_terms
# gene ID type -- ENSEMBL, SYMBOL, ENTREZ, etc
gt <- opt$geneIDtype
# User-specified DE gene p-value cutoff (default = 0.05)
pv <- opt$FDRGenes
# Ontology term p-value cutoff (default = 0.05)
pv_ont <- opt$pvalueOnt
# User-specified DE gene log2 fold change cutoff (default = 1X)
fold_change <-opt$foldchange
# Number of DEGs to use for ontology analysis (default = top 200 genes) 
num_genes <- opt$numDEGs
# Gene Ontology branch - BP, CC, MF or ALL (default = BP)
go_ont <- opt$gene_ontology
# Species
species = opt$species
species = gsub("_"," ",species)

print("DEG file")
DEfile
print("Number of Ontology terms to report")
num_terms
print("gene ID type -- ENSEMBL, SYMBOL, ENTREZ, etc")
gt
print("User-specified DE gene p-value cutoff (default = 0.05)")
pv
print("Ontology term p-value cutoff (default = 0.05)")
pv_ont
print("User-specified DE gene log2 fold change cutoff (default = 1X)")
fold_change
print("Number of DEGs to use for ontology analysis (default = top 200 genes)")
num_genes
print("Gene Ontology branch - BP, CC, MF or ALL (default = BP)")
go_ont
print("Species")
species


#------------species in molecular signatures databases-----------#
dbs <- read.delim("species_databases.txt",sep="\t",row.names = "msig_species")
# shows all species
print(rownames(dbs))
# save msig species library as "orgdb"
orgdb <- as.character(dbs[rownames(dbs)==species,]$db)
# save kegg species identifier as "kegg_sp" (used in )
kegg_sp <- as.character(dbs[rownames(dbs)==species,]$KEGG_species)
# save reactome species identifier as "reactome_org"
reactome_org <- as.character(dbs[rownames(dbs)==species,]$reactome_species)

# Install and load the Entrez species genome annotation library for the species; or load it if it is already installed
if (orgdb %in% rownames(installed.packages())){
  library(orgdb,character.only = TRUE)
} else {
  BiocManager::install(orgdb,character.only = TRUE)
  library(orgdb,character.only = TRUE)
}

# Load the edgeR DE results
d1 <- read.delim(DEfile,header = T, sep = '\t', as.is = T)
colnames(d1)[1] <- "gene"
print("The DEG results table you supplied looks like this:")
head(d1)
print("Number of genes in results table: ",str(dim(d1)[1]))
# subset genes on FDR
d1 <- d1[(d1$FDR <= pv),] 

print("Number of gene IDs that fail to convert to Entrez IDs:")
# Add Entrez IDs
d1_ids <- bitr(d1$gene, fromType = gt, toType=c('ENTREZID'), OrgDb=orgdb)
d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = gt)
names(d1_merged)[2] <- 'log2FC'
print("Number of DEGs: ",str(dim(d1_merged)[1]))



# Set R working directory
project_dir <- opt$dir
#setwd(project_dir)
curr_dir <- getwd()
# Set up output directory and filename prefix
out_dir <- glue("{project_dir}/cluster_profiler_results")
# Create output directory
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
# Create KEGG pathview directory in output directory
kegg_path_dir <- glue("{project_dir}/cluster_profiler_results/kegg_pathview/")
if (!dir.exists(kegg_path_dir)) {dir.create(kegg_path_dir, recursive = TRUE)}


# ---------------------------------------------------------------------
# Define image resize function - dotplots
# ---------------------------------------------------------------------

adjust_image_dims_dp2 <- function(enrichment_dp,num_terms_user){
  test <- enrichment_dp$data
    mlendes <- max(str_length(test$Description))
    # calculate width
    new_width=0.09*(mlendes)+2.75
    if (new_width < 7){
      new_width=7}
    # calculate height
    new_height=0.124*(num_terms_user)+2.98
    #new_height=0.124*(num_terms)+2.98
    if (new_height < 3){
      new_height=3}
    x <- list("new_width"=new_width,"new_height"=new_height)
    return(x)
}

adjust_image_dims_emap <- function(enrichment_table){
  test <- enrichment_table@result
  test <- test[test$pvalue <= pv_ont,]
  if (dim(test)[1] >= 1){
    num_terms <- dim(test)[1]
    # calculate width
    new_width=0.102*(num_terms)+6.93
    new_height=0.102*(num_terms)+6.93
    x2 <- list("new_width"=new_width,"new_height"=new_height)
    return(x2)
  }
}



# ---------------------------------------------------------------------
# Run KEGG
# ---------------------------------------------------------------------



fcs = d1_merged$log2FC
# Add the ENTREZ gene IDs onto the fold changes
names(fcs) <- d1_merged$ENTREZID
# Reorder the genes from largest FC to smallest FC
geneList <- sort(fcs,decreasing = T)
# Number of DE genes with FDR  < 0.05:
length(geneList)
# Take DE genes with > 1 fold change
genes <-  names(geneList)[abs(geneList) > fold_change]
# Perform GO, KEGG, and Reactome on top X genes
genes <- genes[1:num_genes]

if (kegg_sp != ''){
kegg_enrich <- enrichKEGG(gene = genes, # filtered DE gene list
                 organism     = kegg_sp,
                 pvalueCutoff = pv_ont,
                 pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                 minGSSize = 10, # min number of genes annotated by ontology in gene set
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
                 qvalueCutoff = 0.2)

print("Number of Enriched KEGG terms:")
print(dim(kegg_enrich)[1])

if (dim(kegg_enrich)[1] >= 1){
 write.table(kegg_enrich,paste0(out_dir,"/Table_KEGG_enrichment_results.txt"),sep="\t",quote=FALSE)

 if (num_terms == Inf){
 nt = dim(kegg_enrich)[1]
 lab="allterms"
 }
 else{
 if (dim(kegg_enrich)[1] <= num_terms){
 nt = dim(kegg_enrich)[1]
 lab=paste0(str(nt),"terms")
 }
 if (dim(kegg_enrich)[1] > num_terms){
 nt = num_terms
 lab=paste0(str(nt),"terms")
 }
 }
 ##----- write out KEGG pathview results----------#
 print("Running KEGG pathview")
 setwd(kegg_path_dir)
 for (kegg_term in kegg_enrich@result$ID[1:nt]){
 kegg_desc <- gsub(" ","_",kegg_enrich@result[kegg_enrich@result$ID==kegg_term,]$Description)
 kegg_desc <- gsub("/","-",kegg_desc)
 pathview(gene.data  = geneList,
 pathway.id = kegg_term,
 species    = kegg_sp,
 out.suffix = kegg_desc,
 limit      = list(gene=max(abs(geneList)), cpd=1))
 }
 setwd(curr_dir)
 ##--------- end of write out pathview results-----#
 dp <- dotplot(kegg_enrich, showCategory=nt) + ggtitle("KEGG Term over representation")
 print(nt)
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_KEGG_over_representation_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 x2 <- adjust_image_dims_emap(kegg_enrich)
 if (dim(kegg_enrich)[1] >= 2){
 pdf(paste0(out_dir,"/emapplot_KEGG_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(emapplot(kegg_enrich,showCategory=nt))
 dev.off()
 }
 kegg_enrich <- setReadable(kegg_enrich,orgdb,'ENTREZID')
 pdf(paste0(out_dir,"/cnetplot_KEGG_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(cnetplot(kegg_enrich,foldChange = geneList,showCategory=nt))
 dev.off()
}
}


# ---------------------------------------------------------------------
# Run GO
# ---------------------------------------------------------------------


if (orgdb != ''){
ego <- enrichGO(gene          = genes,
                #universe      = names(geneList),
                OrgDb         = orgdb,
                ont           = go_ont,
                pAdjustMethod = "BH",
                pvalueCutoff  = pv_ont,
                qvalueCutoff  = 0.05,
		
        	readable      = TRUE)

print("Number of Enriched GO terms:")
print(dim(ego)[1])

if (dim(ego)[1] >= 1){
 write.table(ego,paste0(out_dir,"/Table_GO_enrichment_results.txt"),sep="\t",quote=FALSE)
 if (num_terms == Inf){
 nt = dim(ego)[1]
 lab="allterms"
 }
 else{
 if (dim(ego)[1] <= num_terms){
 nt = dim(ego)[1]
 lab=paste0(str(nt),"terms")
 }
 if (dim(ego)[1] > num_terms){
 nt = num_terms
 lab=paste0(str(nt),"terms")
 }
 }
 print(nt)
 print(lab)
 dp <- dotplot(ego, showCategory=nt) + ggtitle("GO Term over representation")
 #x <- adjust_image_dims_dp(ego,nt)
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_GO_over_representation_",go_ont,"_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 x2 <- adjust_image_dims_emap(ego)
 if (dim(ego)[1] >= 2){
 pdf(paste0(out_dir,"/emapplot_GO_",go_ont,"_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(emapplot(ego,showCategory=nt))
 dev.off()
 }
 pdf(paste0(out_dir,"/cnetplot_GO_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(cnetplot(ego,foldChange = geneList,showCategory=nt))
 dev.off()
}
}



# ---------------------------------------------------------------------
# Run Reactome
# ---------------------------------------------------------------------

if (reactome_org != ''){
reactome_enrich <- enrichPathway(gene = genes, # filtered DE gene list
                 organism     = reactome_org,
                 pvalueCutoff = pv_ont,
                 pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                 minGSSize = 5, # min number of genes annotated by ontology in gene set
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
                 qvalueCutoff = 0.2,
		 readable = TRUE)

print("Number of Enriched Reactome terms:")
print(dim(reactome_enrich)[1])

if (dim(reactome_enrich)[1] >= 1){
 write.table(reactome_enrich,paste0(out_dir,"/Table_Reactome_enrichment_results.txt"),sep="\t",quote=FALSE)
 if (num_terms == Inf){
 nt = dim(reactome_enrich)[1]
 lab="allterms"
 }
 else{
 nt = num_terms
 lab=paste0(str(nt),"terms")
 }
 dp <- dotplot(reactome_enrich, showCategory=nt) + ggtitle("Reactome Term over representation")
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_Reactome_over_representation_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 x2 <- adjust_image_dims_emap(reactome_enrich)
 if (dim(reactome_enrich)[1] >= 2){
 pdf(paste0(out_dir,"/emapplot_Reactome_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(emapplot(reactome_enrich,showCategory=nt))
 dev.off()
 }
 pdf(paste0(out_dir,"/cnetplot_Reactome_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(cnetplot(reactome_enrich,foldChange = geneList,showCategory=nt))
 dev.off()
}
}



setwd(project_dir)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")




