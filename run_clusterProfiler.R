#!/usr/bin/env Rscript

# Updates to clusterProfiler wrapper:
# v1.3
# Added option to enable/disable kegg pathview (default = FALSE) since it takes so long to run
# Writes out parameter file with all the clusterProfiler wrapper options that were used
# cnet and emap plot sizes scale better with number of terms
 
# v1.4
# runs GSEA (default GSEA disabled, default MsigDB = hallmark database)

# v1.5
# GOplot is available for GO ORA
# New package additions: "tibble","msigdbr","openxlsx"

# v1.6
# GSEA preranking method option is added



#######################################################################
# Install and load required R packages
#######################################################################

req_pkgs <- list("tidyverse","clusterProfiler","RColorBrewer","glue","msigdbr","pathview","readr","stringr","enrichplot","ggplot2","ReactomePA","tibble","msigdbr","openxlsx","DOSE")

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
              help="Biological process (BP), Cellular Component (CC), Molecular function (MF), or all ontologies (ALL)",metavar="character"),
 make_option("--keggpathview",
	      default=FALSE,
	      help="Enable KEGG pathview plots for KEGG terms (TRUE/FALSE)"),
 make_option("--gsea",
	      default=FALSE,
	      help="Enable GSEA using msigdb (TRUE/FALSE). Default Molecular Signatures Database (msigdb) for GSEA = Hallmark (h). See --msigDBs option for other available databases."),
 make_option("--msigdbs",
	      type="character",
	      default="h",
	      help="Molecular Signatures Database (msigdb) for GSEA or comma-separated list of databases or use 'all' to run through all available databases for the species. Otherwise default database is hallmark (h). Available databases: c1, c2_cgp, c2_cp, c2_cp_biocarta, c2_cp_kegg, c2_cp_pid, c2_cp_reactome, c2_cp_wikipathways, c3_mir_mir_legacy, c3_mir_mirdb,c3_tft_gtrd, c3_tft_tft_legacy, c4_cgn, c4_cm, c5_go_bp, c5_go_cc, c5_go_mf, c5_hpo, c6, c7, c8, h. For more information, see: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp"),
 make_option("--prerankMethod",
              type="character",
              default="pval",
	      help="Method of ranking gene list for GSEA. pval will sort by -log10(pval) * (direction of fold change), and lfc will sort by log2(fold change). Note: if FDR column is not provided, it will automatically sort by log2FC. Check the parameter_setting.txt file to check the method run.")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("DEG results table file is not specified (-f).", call.=FALSE)
}

if (is.null(opt$species)){
  print_help(opt_parser)
  stop("No species is specified (-s). Provide genus species name with underscores. eg. Mus_musculus or Canus_lupus_familiaris.", call.=FALSE)
}

#######################################################################
# Set parameters
#######################################################################

# DEG results file
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
# Enable KEGG Pathview
enable_KEGG_pathview <- opt$keggpathview
# Enable GSEA
enable_GSEA <- opt$gsea
# MsigDBs
msigdbs_string <- opt$msigdbs
if (grepl(",",msigdbs_string, fixed = TRUE)){
msigdbs_list <- str_split(msigdbs_string,',')[[1]]
}else{
msigdbs_list <- msigdbs_string
}
print(enable_GSEA)
print(msigdbs_string)
input_gsea_sorting_method=opt$prerankMethod


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
print("enable_KEGG_pathview")
enable_KEGG_pathview
print("enable_GSEA")
enable_GSEA
print("msigdbs_string")
msigdbs_string


# save parameters
parameters_towrite <- data.frame(parameter=as.character(),setting=as.character())
parameters_towrite <- rbind(parameters_towrite,c("DEG file:",DEfile))
parameters_towrite <- rbind(parameters_towrite,c("Species:",species))
parameters_towrite <- rbind(parameters_towrite,c("User-specified DE gene p-value cutoff (default = 0.05):",pv))
parameters_towrite <- rbind(parameters_towrite,c("User-specified DE gene log2 fold change cutoff (default = 1X):",fold_change))
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs to use for ontology analysis after applying FDR and log2FC cutoffs (default = top 200 genes):",num_genes))
parameters_towrite <- rbind(parameters_towrite,c("Gene ID type (SYMBOL, ENSEMBL, ENTREZ, etc) in DEG file:",gt))
parameters_towrite <- rbind(parameters_towrite,c("Number of Ontology terms to report:",num_terms))
parameters_towrite <- rbind(parameters_towrite,c("Gene Ontology branch to use. Either BP, CC, MF or ALL (default = BP):",go_ont))
parameters_towrite <- rbind(parameters_towrite,c("KEGG pathview enabled:",enable_KEGG_pathview))
parameters_towrite <- rbind(parameters_towrite,c("GSEA enabled:",enable_GSEA))
parameters_towrite <- rbind(parameters_towrite,c("MsigDB for GSEA:",msigdbs_string))
colnames(parameters_towrite) <- c("parameter","setting")

#######################################################################
# Species - set up species databases for GO, KEGG, Reactome, msigdb and load organism R library
#######################################################################

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

#######################################################################
# Load DEG results table
#######################################################################
# Load the edgeR DEG results
d1 <- read.delim(DEfile,header = T, sep = '\t', as.is = T)
colnames(d1)[1] <- "gene"
print("The DEG results table you supplied looks like this:")
head(d1)
print("Number of genes in results table: ",str(dim(d1)[1]))
# subset genes on FDR
d1 <- d1[(d1$FDR <= pv),] 
#print(c("Number of genes after applying ",pv," P-value cutoff, and ",fold_change," log2(Fold change cutoff : ",str(dim(d1)[1])))

print("Number of gene IDs that fail to convert to Entrez IDs:")
# Add Entrez IDs
d1_ids <- bitr(d1$gene, fromType = gt, toType=c('ENTREZID'), OrgDb=orgdb)
d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = gt)
names(d1_merged)[2] <- 'log2FC'
#d1_merged <- d1_merged[(d1_merged$PValue <= opt$pvalueGenes) & (abs(d1_merged$log2FC) > opt$foldchange),]
print("Number of DEGs: ",str(dim(d1_merged)[1]))



#######################################################################
# Set working directory for the analysis
#######################################################################

# Set R working directory
project_dir <- opt$dir
curr_dir <- getwd()
# Set up output directory and filename prefix
out_dir <- glue("{project_dir}/clusterProfiler_results")
# Create output directory
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
# Create KEGG pathview directory in output directory
kegg_path_dir <- glue("{out_dir}/kegg_pathview/")
if (!dir.exists(kegg_path_dir)) {dir.create(kegg_path_dir, recursive = TRUE)}


#######################################################################
# Define dotplot image resizing function
#######################################################################

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


adjust_image_dims_dp <- function(enrichment_table,num_terms_user){
  test <- enrichment_table@result
  test <- test[test$pvalue <= pv_ont,]
  if (dim(test)[1] >= 1){
    if (num_terms_user == Inf){
    num_terms <- dim(test)[1]
    }
    if (num_terms_user > dim(test)[1]){
    num_terms <- dim(test)[1]
    }
    if (num_terms_user < dim(test)[1]){
    num_terms <- num_terms_user 
    }
    mlendes <- max(str_length(test$Description))
    # calculate width
    new_width=0.09*(mlendes)+2.75
    if (new_width < 7){
      new_width=7}
    # calculate height
    new_height=0.124*(num_terms)+2.98
    if (new_height < 3){
      new_height=3}
    x <- list("new_width"=new_width,"new_height"=new_height)
    return(x) 
  }
}


adjust_image_dims_emap <- function(enrichment_table){
  test <- enrichment_table@result
  test <- test[test$pvalue <= pv_ont,]
  if (dim(test)[1] >= 1){
    num_terms <- dim(test)[1]
    # calculate width
    new_width=0.05*(num_terms)+6.93 # old scaling factor = 0.105
    new_height=0.05*(num_terms)+6.93 # old scaling factor = 0.105
    x2 <- list("new_width"=new_width,"new_height"=new_height)
    return(x2)
  }
}


#######################################################################
# Analysis
#######################################################################

# ---------------------------------------------------------------------
# Prep data for GSEA
# ---------------------------------------------------------------------

# define column swap function
swap <- function(DF, n, m)
{
  return (DF[ if (n==m) 1:length(DF) else c( (if (min(n,m)==1) c() else 1:(min(n,m)-1) ), (if (min(n,m)+1 == max(n,m)) (min(n,m)+1):(max(n,m)-1) else c( max(n,m), (min(n,m)+1):(max(n,m)-1), min(n,m))), (if (max(n,m)==length(DF)) c() else (max(n,m)+1):length(DF) ) ) ])
}

d1 <- read.delim(DEfile,header = T, sep = '\t', as.is = T)
colnames(d1)[1] <- "gene"
d1_ids <- bitr(d1$gene, fromType = gt, toType=c('ENTREZID'), OrgDb=orgdb)
d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = gt)
names(d1_merged)[2] <- 'log2FC'
print("Number of DEGs: ",str(dim(d1_merged)[1]))
# swap columns to ENTREZID
d1_merged  <- swap(d1_merged,1,dim(d1_merged)[2])
t_d1 <- tibble(d1_merged)

head(t_d1)

# ---------------------------------------------------------------------
# Create pre-ranked gene list for GSEA function
# ---------------------------------------------------------------------

tk_create_prerank <- function(tibble, method = NULL, logFC_col = NULL, padj_col = NULL, sort = TRUE) {
	# method must be one of: "pval_sign" or "lfc"
	# Check input parameters
	if (is.null(logFC_col)) {
		print("You need to provide the column index for log fold change values")
		stop()
	}
	if (is.null(padj_col)) {
		print("You need to provide the column index for adjusted p-values")
		stop()
	}
	# Create prerank table
	if (method == "pval_sign") {
		tibble <- add_column(tibble, "log2FoldChange_sign" = sign(tibble[[logFC_col]]))
		tibble <- add_column(tibble, "minuslog10pval" = -log10(tibble[[padj_col]]))
		# If any are Inf, set to very high value (e.g. 300)
		tibble$minuslog10pval[which(tibble$minuslog10pval == "Inf")] <- 300
		# Create the prerank_metric: -log10pval * fold change sign
		tibble <- add_column(tibble, "prerank_metric" = (tibble$minuslog10pval * tibble$log2FoldChange_sign))
		# Remove genes that have p values with "NA"
		tibble <- tibble[!(is.na(tibble$minuslog10pval)), ]
		if (sort == TRUE) {
			tibble <- tibble[order(tibble$prerank_metric, decreasing = TRUE), ]
		}
		return(tibble)
	} else if (method == "lfc") {
		tibble <- tibble %>% 
					add_column(prerank_metric = tibble[[logFC_col]]) %>%
					drop_na(padj_col) %>%
					filter(ENTREZID != "unknown")
		if (sort == TRUE) {
			tibble <- tibble %>% 
					arrange(desc(tibble[[logFC_col]]))
		}	
		return(tibble)
	} else {
		print("You need to provide a method argument that is either pval_sign or lfc")
		stop()
	}
}



# ---------------------------------------------------------------------
# GSEA dotplot - code hacked from Enrichplot (DOSE)
# ---------------------------------------------------------------------

ep_str_wrap <- function(string, width) {
    x <- gregexpr(' ', string)
    vapply(seq_along(x),
           FUN = function(i) {
               y <- x[[i]]
               n <- nchar(string[i])
               len <- (c(y,n) - c(0, y)) ## length + 1
               idx <- len > width
               j <- which(!idx)
               if (length(j) && max(j) == length(len)) {
                   j <- j[-length(j)]
               }
               if (length(j)) {
                   idx[j] <- len[j] + len[j+1] > width
               }
               idx <- idx[-length(idx)] ## length - 1
               start <- c(1, y[idx] + 1)
               end <- c(y[idx] - 1, n)
               words <- substring(string[i], start, end)
               paste0(words, collapse="\n")
           },
           FUN.VALUE = character(1)
    )
}

default_labeller <- function(n) {
    function(str){
        str <- gsub("_", " ", str)
        ep_str_wrap(str, n)
    }
}



GSEA_dotplot <- function(gsea_results_df,num_cat2){
label_format=30
label_func <- default_labeller(label_format)
if(is.function(label_format)){
    label_func <- label_format}
gene_lists <- str_split(gsea_results_df$core_enrichment,"/")
num_core_enriched <- c()
for (i in 1:length(gene_lists)){
        num_core_enriched <- c(num_core_enriched, length(gene_lists[[i]]))
}
gsea_results_df$num_core_enriched <- num_core_enriched
gsea_results_df$GeneRatio <- as.numeric(gsea_results_df$num_core_enriched/gsea_results_df$setSize)
gsea_results_df$sign <- "suppressed"
if ("TRUE" %in% (gsea_results_df$NES > 0)){
gsea_results_df[gsea_results_df$NES > 0,]$sign <- "activated"
}
gsea_results_df$Description <- as.character(gsea_results_df$ID)
print(gsea_results_df)
font.size=12
color="p.adjust"
colorBy <- "p.adjust"
idx <- order(gsea_results_df$GeneRatio, decreasing = T)
#idx <- order(gsea_results_df[gsea_results_df$GeneRatio], decreasing = F)
gsea_results_df$Description <- factor(gsea_results_df$Description,levels=rev(unique(gsea_results_df$Description[idx])))
gplot <- ggplot(gsea_results_df, aes_string(x="GeneRatio", y="Description", size="num_core_enriched", color=colorBy)) +
geom_point() + scale_color_continuous(low="red", high="blue", name = color,guide=guide_colorbar(reverse=TRUE)) +
scale_y_discrete(labels = label_func) + ylab(NULL) + theme_dose(font.size) + scale_size(range=c(3, 8))+ facet_wrap(~sign)
return(gplot)
}




# ---------------------------------------------------------------------
# GSEA function
# ---------------------------------------------------------------------


tk_enricher_gsea <- function(sorted_prerank_df, msigdb_long, qvalue_threshold = 0.05, database_types = levels(fmsigdb_long$database), num_cat) {

        prerank_gene_list <- sorted_prerank_df$prerank_metric
        names(prerank_gene_list) <- as.character(sorted_prerank_df$ENTREZID)

        # DATABASES
        # all:
        # database_types <- levels(msigdb_long$database)
        # OR MANUAL EDIT TO CHOOSE FOR YOUR FAVORITE DATABASES
        # database_types <- levels(msigdb_long$database)[grepl("_all", levels(msigdb_long$database))]
        # database_types <- c("c2_cgp", "c2_cp_reactome", "c2_cp_biocarta", "c2_cp_kegg", "h", "c6")


        gsea_msigdb_df <- data.frame()
        gsea_msigdb_list <- list()
        for (dbt in seq_along(database_types)) {
                # get only one database type at a time
                # Print database name
                print(database_types[dbt])
                msigdb_long_curr <- filter(msigdb_long, database == database_types[dbt])
                # Print the number of gene sets in database
                print(length(levels(factor(msigdb_long_curr$ont))))
		set.seed(123)
                gsea_msigdb_curr <- clusterProfiler::GSEA(prerank_gene_list, TERM2GENE = msigdb_long_curr, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.5, seed=TRUE, pAdjustMethod = "BH")
                if (length(gsea_msigdb_curr$ID) > 0) {
                        gsea_msigdb_curr <- setReadable(gsea_msigdb_curr, OrgDb = orgdb, keyType = "ENTREZID")
                        gsea_msigdb_curr@result$gene_set_database <- database_types[dbt]
                        gsea_msigdb_curr@result$gene_set_description <- msigdb_long$description[match(gsea_msigdb_curr$ID, msigdb_long$ont)]
                        # Sometimes with only 1 signif geneset, the qvalues are logical with NA. This causes problems when plotting
                        if (!is.logical(gsea_msigdb_curr@result$qvalues)) {
                                gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$qvalues <= qvalue_threshold, ]
                                gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$qvalues, decreasing = FALSE), ]
                        } else {
                                # Else, filter on adjusted p values (which are similar)
                                gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$p.adjust <= qvalue_threshold, ]
                                gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$p.adjust, decreasing = FALSE), ]
                        }
                        # Combine results from curr df into full df
                        gsea_msigdb_df <- rbind(gsea_msigdb_df, gsea_msigdb_curr_signif_df)
                        # Combine results from curr database to all others
                        gsea_msigdb_list[[database_types[dbt]]] <- gsea_msigdb_curr


                        ##------------save dotplot of signif gsea results and image height and width info------------------##
                        if (dim(gsea_msigdb_curr_signif_df)[1] >= 1){
                        if (num_cat == Inf){
                        nt = length(gsea_msigdb_curr_signif_df$ID)
                        print(nt)
                        }else{if (length(gsea_msigdb_curr_signif_df$ID) >= num_cat){
                        nt = num_cat
			gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[1:num_cat,]}
                        if (length(gsea_msigdb_curr_signif_df$ID) < num_cat) {
                        nt = length(gsea_msigdb_curr_signif_df$ID)}
                        }
			lenchar <- max(str_length(as.character(gsea_msigdb_curr_signif_df$Description)))
			new_w <- (0.196*(lenchar))+1.569
			if(new_w < 10) {new_w <- 10}
			new_h <- 0.24*(nt)+5.04
			print(new_h)
			print(new_w)
			if (!dir.exists(paste0("clusterProfiler_results/msigdb_gsea/", database_types[dbt]))) {dir.create(paste0("clusterProfiler_results/msigdb_gsea/", database_types[dbt]), recursive = TRUE)}
			pdf(paste0("clusterProfiler_results/msigdb_gsea/Dotplot_GSEA_",database_types[dbt],".pdf"), width = new_w, height = new_h)
			print(GSEA_dotplot(gsea_msigdb_curr_signif_df,nt))
			Sys.sleep(10)
			dev.off()
                        ##-------------------------------------------------------------------##
                        }
                        # Make individual plots for all signif gene sets
                        for (j in seq_along(gsea_msigdb_curr_signif_df$ID)) {
                                if (gsea_msigdb_curr_signif_df$NES[j] > 0) {
                                        nes_direction <- "pos"
                                } else {
                                        nes_direction <- "neg"
                                }
                                #if (!dir.exists(paste0("clusterProfiler_results/msigdb_gsea/", database_types[dbt]))) {dir.create(paste0("clusterProfiler_results/msigdb_gsea/", database_types[dbt]), recursive = TRUE)}
                                pdf(paste0("clusterProfiler_results/msigdb_gsea/", database_types[dbt], "/", str_pad(j, 4, pad = "0"), "_", nes_direction, "_", gsea_msigdb_curr_signif_df$ID[j], ".pdf"))
                                print(gseaplot(gsea_msigdb_curr, gsea_msigdb_curr_signif_df$ID[j], title = paste0(gsea_msigdb_curr_signif_df$ID[j], "\n", "NES = ", round(gsea_msigdb_curr_signif_df$NES[j], 4), "\n", "qvalue = ", round(gsea_msigdb_curr_signif_df$qvalues[j], 4))))
                                dev.off()
                        }
			
			#pdf(paste0("clusterProfiler_results/msigdb_gsea/Dotplot_GSEA_",database_types[dbt],".pdf"), width = new_width, height = new_height)
			#print(ggplot(gsea_msigdb_curr_signif_df, aes_string(x="GeneRatio", y="Description", size="num_core_enriched", color=colorBy)) +
                        #geom_point() + scale_color_continuous(low="red", high="blue", name = color,guide=guide_colorbar(reverse=TRUE)) +
                        #scale_y_discrete(labels = label_func) + ylab(NULL) + theme_dose(font.size) + scale_size(range=c(3, 8))+ facet_wrap(~sign))
			#dev.off()

                }
        }
        # Clean up and sort data frame, export to Excel file
        gsea_msigdb_df <- gsea_msigdb_df[order(gsea_msigdb_df$qvalues, decreasing = FALSE), ]
        gsea_msigdb_df$Description <- NULL
        class(gsea_msigdb_df$gene_set_description) <- "hyperlink"

        prerank_excel <- as.data.frame(sorted_prerank_df)
        #print(head(prerank_excel))
        #print(table(prerank_excel$prerank_metric))
        #class(prerank_excel$ncbi_gene_db) <- "hyperlink"
        #class(prerank_excel$ensembl_gene_db) <- "hyperlink"
        prerank_excel_short <- prerank_excel[, c("gene", "prerank_metric")]


        excel_list <- list(gsea_msigdb_df, prerank_excel_short, prerank_excel)
        names(excel_list) <- c("gsea_signif", "preranked_gene_order", "preranked_all_columns")
        wb <- openxlsx::write.xlsx(excel_list, file = paste0("clusterProfiler_results/msigdb_gsea/msigdb_gsea.xlsx"))
        setColWidths(wb, sheet = 1, cols = 1, widths = 60)
        saveWorkbook(wb, paste0("clusterProfiler_results/msigdb_gsea/msigdb_gsea.xlsx"), overwrite = TRUE)

        # Create dir for text file output
        if (!dir.exists("clusterProfiler_results/suppl_files")) {dir.create("clusterProfiler_results/suppl_files")}
        write_tsv(gsea_msigdb_df, "clusterProfiler_results/suppl_files/gsea_msigdb.tsv")
        write_tsv(prerank_excel_short, "clusterProfiler_results/suppl_files/prerank_less_cols.tsv")
        write_tsv(prerank_excel, "clusterProfiler_results/suppl_files/prerank.tsv")
        return(gsea_msigdb_list)
}





# ---------------------------------------------------------------------
# End of GSEA function
# ---------------------------------------------------------------------


# retrieve all gene sets associated with the species
msigdb <- msigdbr(species = species)
# convert colons to underscores in gene set subcategory column
colon_to_underscore <- gsub(":", "_", msigdb$gs_subcat)
database <- tolower(paste0(msigdb$gs_cat, "_", colon_to_underscore))
database <- sub("_$","",database) # MM modification -- this is needed; otherwise the "hallmark gene set" shows 0 genes when you run the cluster profiler wrapper function
table(database)
msigdb$database <- database
# add links
msigdb$description <- paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", msigdb[[1]])
# Create a data frame for clusterProfiler (for genes as Entrez Gene IDs).
msigdb_long <- msigdb %>%
dplyr::select(gs_name, entrez_gene, human_gene_symbol, database, description) %>%
as.data.frame()
# change column names
colnames(msigdb_long) <- c("ont", "gene", "symbol", "database", "description")
# write TSV file of msigdb
# write_tsv(msigdb_long, paste0(out_dir, "/", prefix, "msigdb_long.tsv"))
msigdb_long$database <- factor(msigdb_long$database)



if (NA %in% match(msigdbs_list,c("c1", "c2_cgp", "c2_cp", "c2_cp_biocarta", "c2_cp_kegg", "c2_cp_pid", "c2_cp_reactome", "c2_cp_wikipathways", "c3_mir_mir_legacy", "c3_mir_mirdb","c3_tft_gtrd", "c3_tft_tft_legacy", "c4_cgn", "c4_cm", "c5_go_bp", "c5_go_cc", "c5_go_mf", "c5_hpo", "c6", "c7", "c8", "h","all"))) {
print("ERROR: MSigDB entry is not valid. Check the spelling of MSig databases.")
quit(status=1)
}

if (msigdbs_list == "all"){
 msigdbs_list <- unique(msigdb_long$database)

}

if (input_gsea_sorting_method == "pval"){
if (c("FDR","log2FC") %in% colnames(t_d1)){
method="pval_sign"
}else{
method="lfc"
}
}else{
method="lfc"
}

parameters_towrite <- rbind(parameters_towrite,c("GSEA ranking method:",method))
sorted_preranked_df <- tk_create_prerank(t_d1,method=method,logFC_col = which(colnames(t_d1) == "log2FC"),padj_col = which(colnames(t_d1) == "FDR"))



print(msigdbs_list)
# Run GSEA
tk_enricher_gsea(sorted_preranked_df, msigdb_long, qvalue_threshold = 0.05, database_types = msigdbs_list, num_cat=num_terms)




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
# determine the number of DEGs used for ontology analyses
if (length(genes) <= num_genes){ 
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",length(genes)))
}else{
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",as.character(num_genes)))
}
# Perform GO, KEGG, and Reactome on top X genes
genes <- genes[1:num_genes]

write.table(parameters_towrite,paste0(out_dir,"/Parameter_settings.txt"),sep="\t",quote=FALSE,row.names=FALSE)

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
 if (enable_KEGG_pathview == TRUE){
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
 }
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
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_GO_over_representation_",go_ont,"_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 # make emapplot
 x2 <- adjust_image_dims_emap(ego)
 if (dim(ego)[1] >= 2){
 pdf(paste0(out_dir,"/emapplot_GO_",go_ont,"_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(emapplot(ego,showCategory=nt))
 dev.off()
 }
 # make cnetplot
 pdf(paste0(out_dir,"/cnetplot_GO_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(cnetplot(ego,foldChange = geneList,showCategory=nt))
 dev.off()
 # make goplot
 pdf(paste0(out_dir,"/goplot_GO_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(goplot(ego,showCategory = nt))
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
 #reactome_enrich <- setReadable(reactome_enrich,orgdb,'ENTREZID')
 pdf(paste0(out_dir,"/cnetplot_Reactome_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
 print(cnetplot(reactome_enrich,foldChange = geneList,showCategory=nt))
 dev.off()
}
}


writeLines(capture.output(sessionInfo()), paste0(out_dir,"/sessionInfo.txt"))






