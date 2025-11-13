#!/usr/bin/env Rscript

# clusterProfiler Wrapper version 2.2

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

# v1.7
# added checks to ensure enrichment objects do not return NULL statements.

# v1.8 
# Added data set screening functions. Add try statements to try GSEA with each msigdb. Adds "ggnewscale" package dependency for enrichplot. 
# Fixed KEGG enrich issues for D. melanogaster

# v2.0.1
# Fixed issues with pairwise_termsim() failure by adding try statements

# v2.1
# Development version of clusterProfiler vXXX. Kegg_enrich functions now work
# Added universe for GO/KEGG
# Use pval instead of FDR as default for GSEA 
#	[] Change in the options menu
#	[] 

# v2.1-rocky - may 3, 2024
# Modified .libPath to new location
# Writes out RDS objects for the KEGG, GO, and Reactome ORAs
# Added a try() statement around the kegg pathview function (https://support.bioconductor.org/p/9154314/)


# v2.2 October 28,2025 
# Modified the R .libPath to new location (/common/software/)
# GSEA prerank changes: ranking by p-value ranks genes first by -log10(Pvalue), but if there are ties in rank, will then sort genes by log2FC magnitude.
# 

#.libPaths(c("/panfs/jay/groups/0/riss/public/clusterProfilerWrapper2.1_rocky8/R/x86_64-pc-linux-gnu-library/4.3",.libPaths()))
.libPaths(c("/common/software/install/migrated/clusterprofiler_wrapper/R/x86_64-pc-linux-gnu-library/4.3/",.libPaths()))

#######################################################################
# Parse options
#######################################################################

# Clear out any objects from workspace
rm(list = ls(all.names = TRUE))


# check if optparse is installed
req_pkgs <- list("optparse")

for (i in req_pkgs){
  if (i %in% rownames(installed.packages())){
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }else{
    BiocManager::install(i,character.only = TRUE)
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }
}



library("optparse",verbose=FALSE,quietly=TRUE)
option_list = list(
	make_option(c("-d","--dir"),
			type="character",
			default="./clusterProfiler_results/",
			help="project directory",metavar="character"),
	make_option(c("-n","--num_top_terms"),
			type="numeric",
			default="Inf",
			help="Report top N ontology terms",metavar="character"),
	make_option(c("-g","--numDEGs"),
			type="numeric",
			default="500",
			help="Use top X DE genes for Ontology Term enrichment Analysis. Genes with p-value less than 0.05 (or user-specified cutoff) are sorted (largest to smallest) by the absolute value of log2FC. Top 500 genes are kept by default.",metavar="character"),
	make_option(c("-p","--FDRGenes"),
			type="numeric",
			default="0.05",
			help="FDR cutoff used for differentially expressed genes",metavar="character"),
	make_option(c("-a","--pvalueOnt"),
			type="numeric",
			default="0.05",
			help="p-value cutoff used for ontology term enrichment test. Default = 0.05",metavar="character"),
	make_option(c("-c","--foldchange"),
			type="numeric",
			default="1",
			help="absolute log2(fold change) cutoff for differentially expressed genes. Default is log2FC = 1.",metavar="character"),
	make_option(c("-b","--gene_ontology"),
			type="character",
			default="BP",
			help="Biological process (BP), Cellular Component (CC), Molecular function (MF), or all ontologies (ALL)",metavar="character"),
	make_option("--keggpathview",
		    	type = "logical",
			default=FALSE,
			help="Enable KEGG pathview plots for KEGG terms (TRUE/FALSE)"),
	make_option("--gsea",
		    	type = "logical",
			default=FALSE,
			help="Enable GSEA using msigdb (TRUE/FALSE). Default Molecular Signatures Database (msigdb) for GSEA = Hallmark (h). See --msigDBs option for other available databases."),
	make_option("--msigdbs",
			type="character",
			default="h",
			help="Molecular Signatures Database (msigdb) for GSEA or comma-separated list of databases or use 'all' to run through all available databases for the species. Otherwise default database is hallmark (h). Available databases: c1, c2_cgp, c2_cp, c2_cp_biocarta, c2_cp_kegg, c2_cp_pid, c2_cp_reactome, c2_cp_wikipathways, c3_mir_mir_legacy, c3_mir_mirdb, c3_tft_gtrd, c3_tft_tft_legacy, c4_cgn, c4_cm, c5_go_bp, c5_go_cc, c5_go_mf, c5_hpo, c6, c7_immunesigdb, c7_vax, c8, h. For more information, see: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp"),
	make_option("--prerankMethod",
			type="character",
			default="pval",
			help="Method of ranking gene list for GSEA. 'pval' will sort by -log10(pval) * (direction of fold change), and 'lfc' will sort by log2(fold change). Note: if P-value column is not provided, it will automatically sort by log2FC. Check the parameter_setting.txt file to check the method run."),
	make_option("--numPermutations",
			type="numeric",
			default="10000",
			help="Number of times ranked genes are permutated for GSEA. At least 10k permutations are recommended. Default = 10,000.",metavar="character"),
	make_option("--customGSEAgeneSet",
		    	type="character",
			metavar="character",
			help="A 2-column tab-delimited file containing a custom gene set to include alongside the gene sets specified by the --msigsdb option. column 1 should contain the name of the custom gene set. The header of this column should be 'database'. column 2 should contain the genes (in gene SYMBOL format) that are within the custom gene set. The header of column 2 should be set to 'gene'") 
)


# Positional arguments Required: argument 1 = file, argument 2 = species name, argument 3 =  gene ID type
opt_parser = OptionParser(usage = "%prog [options] DEG_file.txt species_name(e.g. Homo_sapiens) GeneID_format(e.g. SYMBOL, ENSEMBL, ENTREZ)", option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = 3,print_help_and_exit=TRUE);

#print(opt)


#######################################################################
# Install and load required R packages
#######################################################################


req_pkgs <- list("tidyverse","clusterProfiler","RColorBrewer","glue","msigdbr","pathview","readr","stringr","enrichplot","ggplot2","ReactomePA","tibble","msigdbr","openxlsx","DOSE","ggnewscale","funr","optparse")

for (i in req_pkgs){
  if (i %in% rownames(installed.packages())){
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }else{
    BiocManager::install(i,character.only = TRUE)
    library(i,character.only = TRUE,verbose=FALSE,quietly=TRUE)
  }
}

options(ggrepel.max.overlaps = Inf)


#######################################################################
# Set parameters
#######################################################################

# DEG results file
DEfile <-opt$args[1]
# Species
species = opt$args[2]
species = gsub("_"," ",species)
# gene ID type -- ENSEMBL, SYMBOL, ENTREZ, etc
gt <- opt$args[3]

# Number of Ontology terms to report
num_terms <- opt$options$num_top_terms
# User-specified DE gene p-value cutoff (default = 0.05)
pv <- opt$options$FDRGenes
# Ontology term p-value cutoff (default = 0.05)
pv_ont <- opt$options$pvalueOnt
# User-specified DE gene log2 fold change cutoff (default = 1X)
fold_change <-opt$options$foldchange
# Number of DEGs to use for ontology analysis (default = top 200 genes) 
num_genes <- opt$options$numDEGs
# Gene Ontology branch - BP, CC, MF or ALL (default = BP)
go_ont <- opt$options$gene_ontology
# Custom GSEA gene set
CUSTOM_GS <- opt$options$customGSEAgeneSet


# Enable KEGG Pathview
enable_KEGG_pathview <- opt$options$keggpathview
# Enable GSEA
enable_GSEA <- opt$options$gsea
# MsigDBs
msigdbs_string <- opt$options$msigdbs
if (grepl(",",msigdbs_string, fixed = TRUE) == TRUE){
msigdbs_list <- str_split(msigdbs_string,',')[[1]]
}else{
msigdbs_list <- msigdbs_string
}
print(enable_GSEA)
print(msigdbs_string)
print(msigdbs_list)
input_gsea_sorting_method=opt$options$prerankMethod
nPerm <- opt$options$numPermutations

print(nPerm)

print("DEG file")
DEfile
print("Number of Ontology terms to report")
num_terms
print("gene ID type -- ENSEMBL, SYMBOL, ENTREZ, etc")
gt
print("User-specified DE gene FDR cutoff (default = 0.05)")
pv
print("Ontology term p-value cutoff (default = 0.05)")
pv_ont
print("User-specified DE gene log2 fold change cutoff (default = 1X)")
fold_change
print("Number of DEGs to use for ontology analysis (default = top 500 genes)")
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
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs to use for ontology analysis after applying FDR and log2FC cutoffs (default = top 500 genes):",num_genes))
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

module_dir <- dirname(sys.script())
dbs <- read.delim(paste0(module_dir,"/species_databases.txt"),sep="\t",row.names = "msig_species")
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

source(paste0(module_dir,"/dataset_screen.R"))
#######################################################################
# Load DEG results table
#######################################################################
# Load the edgeR DEG results
d1 <- read.delim(DEfile,header = T, sep = '\t', as.is = T)
print("The DEG results table you supplied looks like this:")
head(d1)

# Identify Gene, P-Value, and FDR columns
col_indices <- identify_columns(d1)
colnames(d1)[col_indices$Gene] <- "gene"
colnames(d1)[col_indices$LFC] <- "log2FC"
colnames(d1)[col_indices$PVal] <- "FDR"

print("Update DEG column headers:")
head(d1)

# Check log2FC column
validate_lfc(d1$log2FC)
#if (validate_lfc(d1$log2FC) == "TRUE"){
#  stop("Fold change column may be incorrect. Exiting.")
#}

# Check PValue/FDR column. If FALSE- good to go
if (validate_pval(d1$FDR) == "TRUE"){
  stop("P-value column may be incorrect. Exiting.")
}

# Convert Gene ID Format to ENTREZ
d1_ids <- bitr(d1$gene, fromType = gt, toType=c('ENTREZID'), OrgDb=orgdb)
d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = gt)
##print(paste0("Number of DEGs: ",dim(d1_merged)[1]))
##print("Number of gene IDs that fail to convert to Entrez IDs:")

# Save data frame with ENTREZ IDs for GSEA
d1_gsea <- d1_merged
print(paste0("Number of genes in results table= ",dim(d1_merged)[1]))
# subset genes on FDR
d1_merged <- d1_merged[(d1_merged$FDR <= pv),] 
print(paste0("Number of genes in results table after applying FDR cut-off= ",dim(d1_merged)[1]))
#print(d1_merged)

#######################################################################
# Set working directory for the analysis
#######################################################################

# Set R working directory
project_dir <- opt$options$dir
# Set up output directory and filename prefix
out_dir <- project_dir
# Create output directory
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
# Create KEGG pathview directory in output directory
if (enable_KEGG_pathview == TRUE){
kegg_path_dir <- glue("{out_dir}/kegg_pathview/")
if (!dir.exists(kegg_path_dir)) {dir.create(kegg_path_dir, recursive = TRUE)}
}



#######################################################################
# Analysis
#######################################################################

source(paste0(module_dir,"/clusterprofiler_functions.R"))
# ---------------------------------------------------------------------
# Run GSEA
# ---------------------------------------------------------------------

check_gsea <- gsea_ok(d1_gsea,species,col_indices)

if (enable_GSEA == TRUE  & check_gsea == TRUE){
	# Message
	print("Running GSEA....")

	# swap columns to ENTREZID
	d1_gsea  <- swap(d1_gsea,1,dim(d1_gsea)[2])
	t_d1 <- tibble(d1_gsea)
	head(t_d1)

	# retrieve all gene sets associated with the species
	msigdb <- msigdbr(species = species)
	# convert colons to underscores in gene set subcategory column
	colon_to_underscore <- gsub(":", "_", msigdb$gs_subcat)
	database <- tolower(paste0(msigdb$gs_cat, "_", colon_to_underscore))
	database <- sub("_$","",database) # MM modification -- this is needed; otherwise the "hallmark gene set" shows 0 genes when you run the cluster profiler wrapper function
	table(database)
	msigdb$database <- database
	# add links - this doesn't work anymore so we are removing it
	# msigdb$description <- paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", msigdb[[1]])
	# Create a data frame for clusterProfiler (for genes as Entrez Gene IDs).
	msigdb_long <- msigdb %>%
		dplyr::select(gs_name, entrez_gene, human_gene_symbol, database) %>%
		as.data.frame()
	# change column names
	colnames(msigdb_long) <- c("ont", "gene", "symbol", "database")
	# write TSV file of msigdb
	# write_tsv(msigdb_long, paste0(out_dir, "/", prefix, "msigdb_long.tsv"))
	msigdb_long$database <- factor(msigdb_long$database)


	check_msigdb_list <- function(msigdbs_list){
		mcheck <- TRUE
		db_names <- c("c1", "c2_cgp", "c2_cp", "c2_cp_biocarta", "c2_cp_kegg", "c2_cp_pid", "c2_cp_reactome", "c2_cp_wikipathways", "c3_mir_mir_legacy", "c3_mir_mirdb","c3_tft_gtrd", "c3_tft_tft_legacy", "c4_cgn", "c4_cm", "c5_go_bp", "c5_go_cc", "c5_go_mf", "c5_hpo", "c6","c7_immunesigdb", "c7_vax", "c8", "h","all")
		for (i in 1:length(msigdbs_list)){
			if (!msigdbs_list[i] %in% db_names){
			mcheck <- FALSE
			}
		}
		return(mcheck)
	}

	if (check_msigdb_list(msigdbs_list) == FALSE){
	print("ERROR: MSigDB entry is not valid. Check the spelling of Molecular Signatures databases and make sure that there are no spaces between database names. Valid entries include: c1, c2_cgp, c2_cp, c2_cp_biocarta, c2_cp_kegg, c2_cp_pid, c2_cp_reactome, c2_cp_wikipathways, c3_mir_mir_legacy, c3_mir_mirdb,c3_tft_gtrd, c3_tft_tft_legacy, c4_cgn, c4_cm, c5_go_bp, c5_go_cc, c5_go_mf, c5_hpo, c6,c7_immunesigdb, c7_vax, c8, h, all." )
	quit(status=1)
	}


	if ("all" %in% msigdbs_list){
		msigdbs_list <- unique(msigdb_long$database)
	}

	if (input_gsea_sorting_method == "pval"){
		if ("PValue" %in% colnames(t_d1) & "log2FC" %in% colnames(t_d1)){
		method="pval_sign"
		}else{
		method="lfc"
		}
		}else{
		method="lfc"
	}

	# Save the ranking method to the parameter file
	parameters_towrite <- rbind(parameters_towrite,c("GSEA ranking method:",method))

	# Check if a custom gene set file was provided - check that it is not NULL and 
	if (!is.null(CUSTOM_GS) && nzchar(CUSTOM_GS)) {
 	 validated_df <- check_custom_GS_file(CUSTOM_GS,orgdb)
  	# You can then proceed using validated_df
	}else{validated_df <- NULL}
	

	sorted_preranked_df <- tk_create_prerank(t_d1,method=method,logFC_col = which(colnames(t_d1) == "log2FC"),padj_col = which(colnames(t_d1) == "PValue"))
	print(msigdbs_list)
	# Run GSEA
	tk_enricher_gsea(sorted_preranked_df, msigdb_long, qvalue_threshold = 0.05, database_types = msigdbs_list, num_cat=num_terms,project_dir=project_dir,nperm=nPerm,custom_GS=validated_df)
}


# ---------------------------------------------------------------------
# Run KEGG
# ---------------------------------------------------------------------

# OLD
#fcs = d1_merged$log2FC
#names(fcs) <- d1_merged$ENTREZID

# Subset the DEG table to genes with absolute value log2 fold change > threshold
d1_merged_sub <- d1_merged[abs(d1_merged$log2FC)>= fold_change,]
fcs <- d1_merged_sub$log2FC
# Add the ENTREZ gene IDs onto the fold changes
names(fcs) <- d1_merged_sub$ENTREZID

# Set universe from all tested genes
universe <- d1_merged$ENTREZID

print("After applying log2FC threshold, we have these many genes:")
print(length(fcs))
#print(fcs)


# Reorder the genes from largest FC to smallest FC
#geneList <- sort(fcs,decreasing = T)
geneList <- sort(abs(fcs),decreasing = T)
# Number of DE genes with FDR  < 0.05:
length(geneList)
print(geneList)
# Take DE genes with > X log2 fold change
# genes <-  names(geneList)[abs(geneList) > fold_change]

# Sort the genes by absolute fold change
#genes <-  geneList[abs(geneList) > fold_change]
#print(genes)
#genes <- genes[sort(abs(genes),decreasing=T)]
genes <- names(geneList)

# determine the number of DEGs used for ontology analyses
if (length(genes) <= num_genes){ 
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",length(genes)))
}else{
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",as.character(num_genes)))
}

# determine the number of DEGs used for ontology analyses
if (length(genes) <= num_genes){ 
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",length(genes)))
}else{
parameters_towrite <- rbind(parameters_towrite,c("Number of DEGs that were used for the Ontology analyses given the FDR, fold change, and topX gene parameter settings:",as.character(num_genes)))
}
# Perform GO, KEGG, and Reactome on top X genes
genes <- genes[1:num_genes]
genes <- na.omit(genes)
print(genes)
print(paste0("After applying clusterProfiler FDR, log2FC, and topX gene filters, ",length(genes)," genes will be used for GO, KEGG, and Reactome database ORA."))
write.table(parameters_towrite,paste0(out_dir,"/Parameter_settings.txt"),sep="\t",quote=FALSE,row.names=FALSE)

if (length(genes) < 30){
	print("There are less than 30 genes after applying clusterProfiler filters. Will not perform Over-representation Analysis. Please re-evaluate your clusterProfiler FDR and log2FC cutoffs. The default log2FC of 1 may be too stringent for your data set.")
	break
 }else{
print("Running Over-Representation Analysis...")
if (kegg_sp != ''){
    if (kegg_sp == "dme"){
      kegg_enrich <- enrichKEGG(gene = genes, # filtered DE gene list
                                organism     = kegg_sp,
                                pvalueCutoff = pv_ont,
                                pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                                minGSSize = 10, # min number of genes annotated by ontology in gene set
                                maxGSSize = 500, # max number of genes annotated by ontology in gene set
                                qvalueCutoff = 0.2,
                                keyType = "ncbi-geneid")
    }else{
    print("Running KEGG Over-Representation Analysis...")
    if (length(universe) > 5000){
    kegg_enrich <- enrichKEGG(gene = genes, # filtered DE gene list
                 universe = universe,
		 organism     = kegg_sp,
                 pvalueCutoff = pv_ont,
                 pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                 minGSSize = 10, # min number of genes annotated by ontology in gene set
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
		 qvalueCutoff = 0.2)
    }else{
    kegg_enrich <- enrichKEGG(gene = genes, # filtered DE gene list
                 organism     = kegg_sp,
                 pvalueCutoff = pv_ont,
                 pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                 minGSSize = 10, # min number of genes annotated by ontology in gene set
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
                 qvalueCutoff = 0.2)
    }
}
print("Number of Enriched KEGG terms:")
print(dim(kegg_enrich)[1])

if (is.null(kegg_enrich) == FALSE){
  if (dim(kegg_enrich)[1] >= 1){
 write.table(kegg_enrich,paste0(out_dir,"/Table_KEGG_enrichment_results.txt"),sep="\t",quote=FALSE)
 saveRDS(kegg_enrich,paste0(out_dir,"/KEGG_enrichment_results.rds"))
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
 	print("Running KEGG pathview...")
	workdir <- getwd() 
	setwd(kegg_path_dir)
 	for (kegg_term in kegg_enrich@result$ID[1:nt]){
 	kegg_desc <- gsub(" ","_",kegg_enrich@result[kegg_enrich@result$ID==kegg_term,]$Description)
 	kegg_desc <- gsub("/","-",kegg_desc)
	try(pathview(gene.data  = geneList,
	pathway.id = kegg_term,
 	species    = kegg_sp,
 	out.suffix = kegg_desc,
 	limit      = list(gene=max(abs(geneList)), cpd=1)))
 	}
 	setwd(workdir)
 	##--------- end of write out pathview results-----#
 }
 print("Making KEGG ORA Dotplot...")
 dp <- dotplot(kegg_enrich, showCategory=nt,label_format = 300) + ggtitle("KEGG Term over representation")
 print(nt)
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_KEGG_over_representation_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 x2 <- adjust_image_dims_emap(kegg_enrich)
 if (dim(kegg_enrich)[1] >= 2){
 print("Making KEGG ORA Emapplot...")
 try(make_emapplot_cnetplot_kegg(enrichment_object=kegg_enrich,out_dir=out_dir,x2=x2,lab=lab,nt=nt,geneList=geneList,orgdb=orgdb))
 }
}
}
}

# ---------------------------------------------------------------------
# Run GO
# ---------------------------------------------------------------------


if (orgdb != ''){
print("Running GO Over-Representation Analysis...")
if (length(universe) > 5000){
ego <- enrichGO(gene          = genes,
                universe      = universe,
                OrgDb         = orgdb,
                ont           = go_ont,
                pAdjustMethod = "BH",
                pvalueCutoff  = pv_ont,
                qvalueCutoff  = 0.2,
        	readable      = TRUE)
}else{
ego <- enrichGO(gene          = genes,
                OrgDb         = orgdb,
                ont           = go_ont,
                pAdjustMethod = "BH",
                pvalueCutoff  = pv_ont,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
}
print("Number of Enriched GO terms:")
print(dim(ego)[1])

if (is.null(ego) == FALSE){
if (dim(ego)[1] >= 1){
 write.table(ego,paste0(out_dir,"/Table_GO_enrichment_results.txt"),sep="\t",quote=FALSE)
 saveRDS(ego,paste0(out_dir,"/GO_enrichment_results.rds"))
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
 print("Making GO ORA Dotplot...")
 dp <- dotplot(ego, showCategory=nt,label_format = 300) + ggtitle(paste0("GO Term over representation - ",length(genes)))
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_GO_over_representation_",go_ont,"_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 # make emapplot,cnetplot,goplot
 x2 <- adjust_image_dims_emap(ego)
 if (dim(ego)[1] >= 2){
 print("Making GO ORA Emapplot...")
 GO_emapplot_TryCatch(enrichment_object=ego,out_dir=out_dir,go_ont=go_ont,x2=x2,lab=lab,nt=nt)
 }
}
}
}



# ---------------------------------------------------------------------
# Run Reactome
# ---------------------------------------------------------------------

if (reactome_org != ''){
print("Running Reactome Over-Representation Analysis...")
if (length(universe) > 5000){
reactome_enrich <- enrichPathway(gene = genes, # filtered DE gene list
                 organism     = reactome_org,
		 universe = universe,
		 pvalueCutoff = pv_ont,
                 pAdjustMethod = "BH", # Benjamini-Hochberg correction of p-value (false discovery rate)
                 minGSSize = 5, # min number of genes annotated by ontology in gene set
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
                 qvalueCutoff = 0.2,
                 maxGSSize = 500, # max number of genes annotated by ontology in gene set
                 qvalueCutoff = 0.2,
                 readable = TRUE)
}

print("Number of Enriched Reactome terms:")
print(dim(reactome_enrich)[1])

if (is.null(reactome_enrich) == FALSE){
 if (dim(reactome_enrich)[1] >= 1){
 write.table(reactome_enrich,paste0(out_dir,"/Table_Reactome_enrichment_results.txt"),sep="\t",quote=FALSE)
 saveRDS(reactome_enrich,paste0(out_dir,"/Reactome_enrichment_results.rds"))
 if (num_terms == Inf){
 nt = dim(reactome_enrich)[1]
 lab="allterms"
 }
 else{
 nt = num_terms
 lab=paste0(str(nt),"terms")
 }
 dp <- dotplot(reactome_enrich, showCategory=nt,label_format = 300) + ggtitle(paste0("Reactome Term over representation - ", length(genes)))
 x <- adjust_image_dims_dp2(dp,nt)
 pdf(paste0(out_dir,"/Dotplot_Reactome_over_representation_",lab,".pdf"), width = x$new_width, height = x$new_height)
 print(dp)
 dev.off()
 x2 <- adjust_image_dims_emap(reactome_enrich)
 if (dim(reactome_enrich)[1] >= 3){
 try(make_emapplot_cnetplot_reactome(enrichment_object=reactome_enrich,out_dir=out_dir,x2=x2,lab=lab,nt=nt))
 }
}
}
}

} # end of ifelse statement of < 30 genes

writeLines(capture.output(sessionInfo()), paste0(out_dir,"/sessionInfo.txt"))




