# clusterProfiler_wrapper 1.0


See the “clusterProfiler_wrapper_intro.pdf" for more details on the clusterProfiler wrapper and for examples on how to run the wrapper.

Important default behavior of the wrapper script:
1. Will subset an edgeR differential expression results table to the genes that have FDR < 0.05 (can adjust this with “-p” option) and abs(log2FC)  > 1 (can adjust this with “-c” option)
2. Will take the top 200 (can adjust this with “-g” option)  genes from (1) to use for GO, KEGG, Reactome term overrepresentation analysis.
3. Then, it will plot all significant terms (p-value < 0.05) for each (can adjust the number of significant terms reported with “-n”).



Feel free to check out clusterprofiler wrapper v1.6 by clicking on the branch link above.
