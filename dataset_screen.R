# Define functions to summarize the dataset and screen it for inappropriate
# inputs to the enrichment and/or GSEA functions.

options(warn=1)

# Identify the columns of interest - gene identifier, adjusted P-value, and
# log(fold change). Takes a data.frame as input and returns a list with three
# elements. If a column is not found, NULL will be returned.
identify_columns <- function(dataset) {
    # A bit silly, but we want to capture many reasonable ways people name the
    # gene ID column. We will generate combinations of "gene" and
    # id/name/symbol with various capitalisations and put together with various
    # delimiters.
    gcol_base <- c("gene", "GENE", "Gene", "genes")
    gcol_suffix <- c("id", "ID", "name", "Name", "NAME", "symbol", "Symbol", "SYMBOL")
    gcol_comb <- expand.grid(gcol_base, gcol_suffix)
    gcol_universe <- c(
        gcol_base,
        gcol_suffix,
        paste(gcol_comb$Var1, gcol_comb$Var2, sep=""),
        paste(gcol_comb$Var1, gcol_comb$Var2, sep="."),
        paste(gcol_comb$Var1, gcol_comb$Var2, sep="_"),
        paste(gcol_comb$Var1, gcol_comb$Var2, sep="-"))
    # Define a similar vector for the P-value
    pval_universe <- c("FDR", "fdr", "padj", "p.adj", "P.Adj", "p.Adj",
                       "p_adj", "P_adj", "P_Adj")
    # And one for the log(fold change)
    lfc_universe <- c("log2FoldChange", "lfc", "LFC", "logFC", "l2FC", "l2fc",
                      "L2FC", "logFoldChange")
    # Check for these names in the colnames of the data.frame
    gcol_in_dat <- which(gcol_universe %in% colnames(dataset))
    pval_in_dat <- which(pval_universe %in% colnames(dataset))
    lfc_in_dat <- which(lfc_universe %in% colnames(dataset))
    # Check that there is one and only one of each of these columns. Print
    # warning messages if not
    if(length(gcol_in_dat) > 1) {
        gcol_use <- gcol_universe[gcol_in_dat[1]]
        warning(paste(
            "More than one gene column found. Using ",
            gcol_use,
            " as the gene identifier column."))
        gcol_idx <- which(colnames(dataset) == gcol_use)
    } else if(length(gcol_in_dat) == 0) {
        gcol_idx <- NULL
        warning("No gene identifier column found!")
    } else {
        gcol_use <- gcol_universe[gcol_in_dat]
        gcol_idx <- which(colnames(dataset) == gcol_use)
    }

    if(length(pval_in_dat) > 1) {
        pval_use <- pval_universe[pval_in_dat[1]]
        warning(paste(
            "More than one P-value column found. Using ",
            pval_use,
            " as the P-value column."))
        pval_idx <- which(colnames(dataset) == pval_use)
    } else if(length(pval_in_dat) == 0) {
        pval_idx <- NULL
        warning("No P-value column found!")
    } else {
        pval_use <- pval_universe[pval_in_dat]
        pval_idx <- which(colnames(dataset) == pval_use)
    }

    if(length(lfc_in_dat) > 1) {
        lfc_use <- lfc_universe[lfc_in_dat[1]]
        warning(paste(
            "More than one log(fold change) found. Using ",
            lfc_use,
            " as the log(fold change) column."))
        lfc_idx <- which(colnames(dataset) == lfc_use)
    } else if(length(lfc_in_dat) == 0) {
        lfc_idx <- NULL
        warning("No log(fold change) column found!")
    } else {
        lfc_use <- lfc_universe[lfc_in_dat]
        lfc_idx <- which(colnames(dataset) == lfc_use)
    }
    # Return the list of column indices
    return(list(Gene=gcol_idx, PVal=pval_idx, LFC=lfc_idx))
}

# Validate the log(fold change) column. It should contain numeric values
# between -Inf and +Inf. Print warnings if there are huge values, though,
# because these are weird. Returns FALSE if there are no issues and returns
# TRUE if there are issues.
validate_lfc <- function(lfc_col) {
    # Check that they are not NA.
    trouble <- FALSE
    if(any(is.na(lfc_col))) {
        trouble <- TRUE
        warning("Some genes have a log(fold change) of NA.")
    }
    # Check that it is numeric
    if(!is.numeric(lfc_col)) {
        trouble <- TRUE
        warning("The log(fold change) column has non-numeric values.")
    }
    # Check for big values
    if(any(lfc_col > 6 | lfc_col < -6)) {
        warning("The log(fold change) column has very large values.")
    }
    return(trouble)
}

# Validate the P-value column. It should contain numeric values between 0 and
# 1.
validate_pval <- function(pval_col) {
    trouble <- FALSE
    if(any(is.na(pval_col))) {
        trouble <- TRUE
        warning("Some genes have a P-value of NA.")
    }
    if(!is.numeric(pval_col)) {
        trouble <- TRUE
        warning("The P-value column has non-numeric values.")
    }
    if(any(pval_col < 0 | pval_col > 1)) {
        trouble <- TRUE
        warning("The P-value column has values outside of [0, 1].")
    }
    return(trouble)
}

# Return a boolean that determines whether we can run GSEA.
gsea_ok <- function(dataset, species, col_idx) {
    ok <- TRUE
    # First, check the number of genes. There should be more than 10000 for
    # those that are not S. cerevisiae, and more than 3000 for yeast.
    if(species == "Saccharomyces_cerevisiae") {
        if(nrow(dataset) < 3000) {
            warning("Less than 3,000 genes in the dataset, GSEA disabled.")
            ok <- FALSE
        }
    } else {
        if(nrow(dataset) < 10000) {
            warning("Less than 10,000 genes in the dataset, GSEA disabled.")
            ok <- FALSE
        }
    }
    # Then, check the presence/absence of the log(fold change) column. If there
    # is no log(fold change) column, then disable GSEA.
    if(is.null(col_idx$LFC)) {
        warning("Could not find log(fold change) column, GSEA disabled.")
        ok <- FALSE
    }
    return(ok)
}
