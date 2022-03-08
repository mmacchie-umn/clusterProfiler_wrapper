# clusterProfiler_wrapper 1.0


#### Required arguments:
1) DEG table in edgeR format (-f). Tab-delimited text file.
2) genus species name (-s). Genus and species should be separated by an underscore ("_")
3) gene name format (-i). Usually in one of the following formats: ENSEMBL, ENTREZ, SYMBOL

### Important. By default, GO, KEGG, and Reactome over-representation analysis (ORA) will perform the following:
1) Will subset the DEG table by FDR (default = 0.05, can set with "-p" option) and fold change (default = 1X, can change with "-c" option)
2) Will use the top 200 genes after sorting genes by absolute fold change (default = 200 genes, can set number of genes with "-g" option) for ORA
3) Will plot ALL significantly enriched ontology terms (can set number of terms to plot with "-n", and it will plot most significant by FDR. Can set the significance threshold for ontology terms with "-a")


# Example commands

### To run top 200 D. melanogaster genes (in gene SYMBOL format) through GO, KEGG, and Reactome over-representation analysis 

Rscript run_clusterProfiler.R -f Drosophila_melanogaster_edgeR_DEG_table.txt -s Drosophila_melanogaster -i SYMBOL 




See the “clusterProfiler_wrapper_intro.pdf" for more details on the clusterProfiler wrapper and for examples on how to run the wrapper.



Feel free to check out clusterprofiler wrapper v1.6 (handles GSEA) by clicking on the branch link above.
