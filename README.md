# ClusterProfiler Wrapper 2.0

Wrapper can now perform GSEA on an edgeR differentially expressed gene results table by including the option: --gsea TRUE. By default, genes will be ranked by their log2FC.
Note: GSEA should only be run on full edgeR results tables that have all tested genes, if less than 10,000 tested genes (or < 3,000 tested genes for S. cerevisiae) are detected in the table, GSEA will be automatically disabled 

#### Required arguments:
1) DEG table in edgeR format (-f). Tab-delimited text file.
2) genus species name (-s). Genus and species should be separated by an underscore ("_")
3) gene name format (-i). Usually in one of the following formats: ENSEMBL, ENTREZ, SYMBOL

### Important. By default, GO, KEGG, and Reactome over-representation analysis (ORA) will perform the following:
1) Will subset the DEG table by FDR (default = 0.05, can set with "-p" option) and fold change (default = 1X, can change with "-c" option)
2) Will use the top 200 genes after sorting genes by absolute fold change (default = 200 genes, can set number of genes with "-g" option) for ORA
3) Will plot ALL significantly enriched ontology terms (can set number of terms to plot with "-n", and it will plot most significant by FDR. Can set the significance threshold for ontology terms with "-a")


# Example commands

### To run top 200 D. melanogaster genes (in gene SYMBOL format) through GO, KEGG, and Reactome over-representation analysis and KEGG pathview analysis

Rscript run_clusterProfiler.R -f Test_Data/Drosophila_melanogaster_edgeR_DEG_table.txt -s Drosophila_melanogaster -i SYMBOL --keggpathview TRUE



# clusterProfiler_wrapper version updates

v2.0
- Added differential gene expression table screening functions
- Added try statements for GSEA function (fixes errors that result from odd DEG tables)
- Fixed enrichplot label issues
- Fixed KEGGenrich gene ID issues for Drosophila melanogaster

v1.6
- GSEA preranking method option is added

v1.5
- GOplot added for GO ORA 

v1.4
- runs GSEA (default GSEA disabled, default MsigDB = hallmark database)

v1.3
- Added option to enable/disable kegg pathview (default = FALSE)
- Writes out parameter file with all the clusterProfiler wrapper options used
- cnet and emap plot sizes scale better with number of terms


