# clusterProfiler_wrapper 1.0


See the “clusterProfiler_wrapper_intro.pdf" for more details on the clusterProfiler wrapper and for examples on how to run the wrapper.
Clusterprofiler is also available 

Important default behavior of the wrapper script:
1. Will subset an edgeR differential expression results table to the genes that have FDR < 0.05 (can adjust this with “-p” option) and log2FC  > 1 (can adjust this with “-c” option)
2. Will take the top 200 (can adjust this with “-g” option)  genes from (1) to use for GO, KEGG, Reactome term overrepresentation analysis.
3. Then, it will plot all significant terms for each (can adjust the number of significant terms reported with “-n”).
