# ClusterProfiler Wrapper 2.0

Wrapper can now perform GSEA on an edgeR differentially expressed gene results table by including the option: --gsea TRUE. By default, genes will be ranked by their log2FC.
Note: GSEA should only be run on full edgeR results tables that have all tested genes, not only the significant ones.


# Example commands



# clusterProfiler_wrapper version updates

v2.0
- Added differential gene expression table screen functions
- Added try statements for GSEA function
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


