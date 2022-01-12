# ClusterProfiler Wrapper 1.6

Wrapper can now perform GSEA on an edgeR differentially expressed gene results table by including the option: --gsea TRUE. By default, genes will be ranked by their log2FC.
Note: GSEA should only be run on full edgeR results tables that have all tested genes, not only the significant ones.


# clusterProfiler_wrapper v1.6 updates

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


