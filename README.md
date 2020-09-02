# graph-generator
Python utility for generating PCA plots and heatmaps from expression data

    usage: python PCA.py [-h] [--exp_set [EXP_SET [EXP_SET ...]]] [--exp_list EXP_LIST]
                  [--metadata_set [METADATA_SET [METADATA_SET ...]]]
                  [--metadata_list METADATA_LIST] [--gene_list GENE_LIST] [-S S]
                  [-T T]


    Example usage: python PCA.py 
    --exp_set GSE125197_rpkm_minimal.txt GSE125197_rpkm_repeated.txt 
    --metadata_set metadata_example.csv 
    -S sample_id -T cell_type 
    --gene_list cell_specific.csv 

Designed for python3. To create a heatmap with a specific set of genes (e.g. cell-type specific genes or any other set of genes), please use the --gene_list parameter.  

Details on input parameters:

        --exp_set [EXP_SET [EXP_SET ...]]
                              Expression matrices, seperated by space. Tab-delimited
                              (.txt) or comma-seperated (.csv) files are accepted
        --exp_list EXP_LIST   Text file of expression matrix list, one file per
                              line. Ignored if --exp_set is also specified
        --metadata_set [METADATA_SET [METADATA_SET ...]]
                              Metadata files with sample id and cell type, seperated
                              by space. Tab-delimited (.txt) or comma-seperated
                              (.csv) files are accepted
        --metadata_list METADATA_LIST
                              Text file of metadata file list, one file per line.
                              Ignored if --metadata_set is also specified
        --gene_list GENE_LIST
                              Cell type specific genes, used for plotting heatmaps.
                              Tab-delimited (.txt) or comma-seperated (.csv), no
                              header. 1st column (required) contains gene
                              identifiers. 2nd column (optional) contains the
                              associated cell type for the gene.
        -S S                  Name of column used to specify the sample name in the
                              metadata files.
        -T T                  Name of column used to specify the cell type in the
                              metadata files.
