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

make_emapplot_cnetplot_kegg <- function(enrichment_object,out_dir,x2,lab,nt,geneList,orgdb){
	enrichment_object <- pairwise_termsim(enrichment_object)
	pdf(paste0(out_dir,"/emapplot_KEGG_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
	print(emapplot(enrichment_object,showCategory=nt))
	dev.off()
        enrichment_object <- setReadable(enrichment_object,orgdb,'ENTREZID')
        pdf(paste0(out_dir,"/cnetplot_KEGG_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
        print(cnetplot(enrichment_object,foldChange = geneList,showCategory=nt))
        dev.off()
	}



make_emapplot_cnetplot_goplot_go <- function(enrichment_object,out_dir,go_ont,x2,lab,nt){
	enrichment_object <- pairwise_termsim(enrichment_object)
	pdf(paste0(out_dir,"/emapplot_GO_",go_ont,"_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
	print(emapplot(enrichment_object,showCategory=nt))
	dev.off()
	pdf(paste0(out_dir,"/cnetplot_GO_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
        print(cnetplot(enrichment_object,foldChange = geneList,showCategory=nt))
        dev.off()
        pdf(paste0(out_dir,"/goplot_GO_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
        print(goplot(enrichment_object,showCategory = nt))
        dev.off()
	}


make_emapplot_cnetplot_reactome <- function(enrichment_object,out_dir,x2,lab,nt){
	enrichment_object <- pairwise_termsim(enrichment_object)
	pdf(paste0(out_dir,"/emapplot_Reactome_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
	print(emapplot(enrichment_object,showCategory=nt))
	dev.off()
        enrichment_object <- setReadable(enrichment_object,orgdb,'ENTREZID')
        pdf(paste0(out_dir,"/cnetplot_Reactome_",lab,".pdf"),height = x2$new_height,width = x2$new_width)
        print(cnetplot(enrichment_object,foldChange = geneList,showCategory=nt))
        dev.off()
	}

GO_emapplot_TryCatch <- function(enrichment_object,out_dir,go_ont,x2,lab,nt){
	tryCatch({ #try this
	make_emapplot_cnetplot_goplot_go(enrichment_object,out_dir,go_ont,x2,lab,nt)
	},
	  # if fails:
 	error=function(error_message){
   	message("There is a problem with pairwise_termsim() at the GO ORA step. Will not export emapplot, cnetplot, and goplot. Here is the R error:")
	message(error_message)
  }
  )
}



#######################################################################
# Define GSEA functions
#######################################################################
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
    tibble <- add_column(tibble, "abs_fc" = abs(tibble[[logFC_col]]))

    # Remove genes that have p values with "NA"
    tibble <- tibble[!(is.na(tibble$minuslog10pval)), ]
    if (sort == TRUE) {
      tibble <- tibble[order(-tibble$prerank_metric, -tibble$abs_fc), ]
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
# Check Custom GSEA file 
# ---------------------------------------------------------------------

check_custom_GS_file <- function(file,orgdb) {
  # Try to read the file, catching errors if file not found
  df <- tryCatch(
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) {
      stop("Custom gene set file not found: ", file)
    }
  )
  # Check number of columns
  if (ncol(df) != 2) {
    stop("GSEA error - the custom gene set file provided must have exactly 2 columns and be tab-delimited. column 1 must have gene set name and column 2 must have the gene IDs in gene SYMBOL format.")
  }

  # Check headers
  if (!all(colnames(df) == c("database", "gene"))) {
    stop("GSEA error - Column headers in the custom gene set file provided. Column 1 header must be 'database' and column 2 header must be 'gene'. They must contain the gene set ID and gene IDs in gene SYMBOL format, respectively.")
  }
  ## Add new columns
  # copy the gene symbol column to "gene" column
  df$symbol <- df$gene
  test <-  bitr(df$symbol, fromType = "SYMBOL", toType=c('ENTREZID'), OrgDb=orgdb)
  df_merge <- merge(df, test, by.x = "gene", by.y = "SYMBOL")
  df <- df_merge[,c(2,4,3)]
  # change the gene set column header to "ont"
  colnames(df)[1:3] <- c("ont","gene","symbol")
  # Create a database column and set the values to "custom"
  df$database <- "custom"
  #df$description <- "custom_gene_set_supplied_by_user"
  ## Add a "Custom_" prefix to the gene set names
  # df$database <- paste0("Custom_",df$database)
  return(df)
}


# ---------------------------------------------------------------------
# GSEA functions
# ---------------------------------------------------------------------

GSEA_TryCatch <- function(prl,mdb,nperm){ 
  tryCatch({
    # try this:
    x <- clusterProfiler::GSEA(prl, TERM2GENE = mdb, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.5, seed=TRUE, pAdjustMethod = "BH",nPermSimple = nperm)
    return(x)
  },
  # if fails:
  error=function(error_message){
    message("And below is the error message from R:")
    message(error_message)
    return(NULL)
  }
  )
}

tk_enricher_gsea <- function(sorted_prerank_df, msigdb_long, qvalue_threshold = 0.05, database_types = levels(msigdb_long$database), num_cat,project_dir,nperm,custom_GS=NULL) {
  
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
#  for (dbt in seq_along(database_types)) {
#    # get only one database type at a time
#    # Print database name
#    print(database_types[dbt])
#    msigdb_long_curr <- filter(msigdb_long, database == database_types[dbt])
#    # Print the number of gene sets in database
#    print(length(levels(factor(msigdb_long_curr$ont))))
  if(is.null(custom_GS)==FALSE){
	msigdb_long_curr <- rbind(custom_GS,msigdb_long)
  }else{msigdb_long_curr <- msigdb_long}

  set.seed(123)
    
    ##gsea_msigdb_curr <- clusterProfiler::GSEA(prerank_gene_list, TERM2GENE = msigdb_long_curr, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.5, seed=TRUE, pAdjustMethod = "BH")
    gsea_msigdb_curr <- GSEA_TryCatch(prerank_gene_list,msigdb_long_curr,nperm)
#    ##if (length(gsea_msigdb_curr$ID) > 0) {
    if (is.null(gsea_msigdb_curr) == FALSE & length(gsea_msigdb_curr$ID) > 0){
      gsea_msigdb_curr <- setReadable(gsea_msigdb_curr, OrgDb = orgdb, keyType = "ENTREZID")
#      gsea_msigdb_curr@result$gene_set_database <- database_types[dbt]
#      gsea_msigdb_curr@result$gene_set_description <- msigdb_long$description[match(gsea_msigdb_curr$ID, msigdb_long$ont)]
      # Sometimes with only 1 signif geneset, the qvalues are logical with NA. This causes problems when plotting
      if (!is.logical(gsea_msigdb_curr@result$qvalue)) {
        gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$qvalue <= qvalue_threshold, ]
        gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$qvalue, decreasing = FALSE), ]
      } else {
        # Else, filter on adjusted p values (which are similar)
        gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$p.adjust <= qvalue_threshold, ]
        gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$p.adjust, decreasing = FALSE), ]
      }
#      # Combine results from curr df into full df
#      gsea_msigdb_df <- rbind(gsea_msigdb_df, gsea_msigdb_curr_signif_df)
#      # Combine results from curr database to all others
#      gsea_msigdb_list[[database_types[dbt]]] <- gsea_msigdb_curr
      
      
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
        if (!dir.exists(paste0(project_dir,"/msigdb_gsea/"))){dir.create(paste0(project_dir,"/msigdb_gsea/"),recursive = TRUE)}
#        if (!dir.exists(paste0(project_dir,"/msigdb_gsea/", database_types[dbt]))) {dir.create(paste0(project_dir,"/msigdb_gsea/", database_types[dbt]), recursive = TRUE)}
#        pdf(paste0(project_dir,"/msigdb_gsea/Dotplot_GSEA_",database_types[dbt],".pdf"), width = new_w, height = new_h)
	pdf(paste0(project_dir,"/msigdb_gsea/Dotplot_GSEA_Results.pdf"), width = new_w, height = new_h)
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
#        #if (!dir.exists(paste0(project_dir,"/msigdb_gsea/", database_types[dbt]))) {dir.create(paste0(project_dir,"/msigdb_gsea/", database_types[dbt]), recursive = TRUE)}
#        pdf(paste0(project_dir,"/msigdb_gsea/", database_types[dbt], "/", str_pad(j, 4, pad = "0"), "_", nes_direction, "_", gsea_msigdb_curr_signif_df$ID[j], ".pdf"))
      	pdf(paste0(project_dir,"/msigdb_gsea/", str_pad(j, 4, pad = "0"), "_", nes_direction, "_", gsea_msigdb_curr_signif_df$ID[j], ".pdf"))
      	print(gseaplot(gsea_msigdb_curr, gsea_msigdb_curr_signif_df$ID[j], title = paste0(gsea_msigdb_curr_signif_df$ID[j], "\n", "NES = ", round(gsea_msigdb_curr_signif_df$NES[j], 4), "\n", "qvalue = ", round(gsea_msigdb_curr_signif_df$qvalue[j], 4))))
        dev.off()
      }
      
    }
    
  if (dim(gsea_msigdb_df)[1] > 0){
  # Clean up and sort data frame, export to Excel file
  gsea_msigdb_df <- gsea_msigdb_df[order(gsea_msigdb_df$qvalue, decreasing = FALSE), ]
#  gsea_msigdb_df$Description <- NULL
#  class(gsea_msigdb_df$gene_set_description) <- "hyperlink"
  
  prerank_excel <- as.data.frame(sorted_prerank_df)
  #print(head(prerank_excel))
  #print(table(prerank_excel$prerank_metric))
  #class(prerank_excel$ncbi_gene_db) <- "hyperlink"
  #class(prerank_excel$ensembl_gene_db) <- "hyperlink"
  prerank_excel_short <- prerank_excel[, c("gene", "prerank_metric")]
  
  
  excel_list <- list(gsea_msigdb_df, prerank_excel_short, prerank_excel)
  names(excel_list) <- c("gsea_signif", "preranked_gene_order", "preranked_all_columns")
  wb <- openxlsx::write.xlsx(excel_list, file = paste0(project_dir,"/msigdb_gsea/msigdb_gsea.xlsx"))
  setColWidths(wb, sheet = 1, cols = 1, widths = 60)
  saveWorkbook(wb, paste0(project_dir,"/msigdb_gsea/msigdb_gsea.xlsx"), overwrite = TRUE)
  
  # Create dir for text file output
  if (!dir.exists(paste0(project_dir,"/suppl_files"))) {dir.create(paste0(project_dir,"/suppl_files"))}
  write_tsv(gsea_msigdb_df, paste0(project_dir,"/suppl_files/gsea_msigdb.tsv"))
  write_tsv(prerank_excel_short, paste0(project_dir,"/suppl_files/prerank_less_cols.tsv"))
  write_tsv(prerank_excel, paste0(project_dir,"/suppl_files/prerank.tsv"))
#  return(gsea_msigdb_list)
  }
}


# define column swap function
swap <- function(DF, n, m)
{
  return (DF[ if (n==m) 1:length(DF) else c( (if (min(n,m)==1) c() else 1:(min(n,m)-1) ), (if (min(n,m)+1 == max(n,m)) (min(n,m)+1):(max(n,m)-1) else c( max(n,m), (min(n,m)+1):(max(n,m)-1), min(n,m))), (if (max(n,m)==length(DF)) c() else (max(n,m)+1):length(DF) ) ) ])
}

