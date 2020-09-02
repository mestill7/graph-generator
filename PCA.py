import pandas as pd
import numpy as np
import sys, getopt
import argparse
# import wget
import os
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns; sns.set()

### Example uses
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='Visualize expression inputs as PCA and heatmaps.',
    epilog='''
    Example usage: python PCA.py 
    --exp_set GSE125197_rpkm_minimal.txt GSE125197_rpkm_repeated.txt 
    --metadata_set metadata_example.csv 
    -S sample_id -T cell_type 
    --gene_list cell_specific.csv \n
    Report issues or feature requests to Github 
    (https://github.com/mestill7/graph-generator)''')
parser.add_argument('--exp_set',nargs='*',type=str, help='Expression matrices, seperated by space. Tab-delimited (.txt) or comma-seperated (.csv) files are accepted')
parser.add_argument('--exp_list',type=str, help='Text file of expression matrix list, one file per line. Ignored if --exp_set is also specified')
parser.add_argument('--metadata_set',nargs='*',type=str, help='Metadata files with sample id and cell type, seperated by space. Tab-delimited (.txt) or comma-seperated (.csv) files are accepted')
parser.add_argument('--metadata_list',type=str, help='Text file of metadata file list, one file per line. Ignored if --metadata_set is also specified')
parser.add_argument('--gene_list',type=str, help='Cell type specific genes, used for plotting heatmaps. Tab-delimited (.txt) or comma-seperated (.csv), no header. 1st column (required) contains gene identifiers. 2nd column (optional) contains the associated cell type for the gene.')
parser.add_argument('-S',type=str, help='Name of column used to specify the sample name in the metadata files.')
parser.add_argument('-T',type=str, help='Name of column used to specify the cell type in the metadata files.')
args = parser.parse_args()
# nargs='*',

# print(args)
arg_dict = vars(args)

if arg_dict['exp_set'] is None and arg_dict['exp_list'] is None:
    print(("No expression matrices provided. Now exiting..."))
    quit()

if arg_dict['metadata_set'] is None and arg_dict['metadata_list'] is None:
    print(("No metadata information provided. Now exiting..."))
    quit()

print(arg_dict)
#Determine file endings for expression set
if arg_dict['exp_set'] is None:
    expression_list = open(arg_dict['exp_list']).read().split("\n") # input file describing the count file names
    expression_list = [x for x in expression_list if x.strip() != '']
    if expression_list[0].split(".")[-1] == "txt":
        rpkm = [pd.read_csv(x,header=0,sep='\t') for x in expression_list]
    if expression_list[0].split(".")[-1] == "csv":
        rpkm = [pd.read_csv(x,header=0,sep=',') for x in expression_list]
else:
    expression_list = [x for x in arg_dict['exp_set'] if x.strip() != '']
    if expression_list[0].split(".")[-1] == "txt":
        rpkm = [pd.read_csv(x,header=0,sep='\t') for x in expression_list]
    if expression_list[0].split(".")[-1] == "csv":
        rpkm = [pd.read_csv(x,header=0,sep=',') for x in expression_list]

if arg_dict['metadata_set'] is None:
    metadata_list = open(arg_dict['metadata_list']).read().split("\n") # input file describing the count file names
    metadata_list = [x for x in metadata_list if x.strip() != '']
    if metadata_list[0].split(".")[-1] == "txt":
        md = [pd.read_csv(x,header=0,sep='\t') for x in metadata_list]
    if expression_list[0].split(".")[-1] == "csv":
        md = [pd.read_csv(x,header=0,sep=',') for x in metadata_list]
else:
    metadata_list = [x for x in arg_dict['metadata_set'] if x.strip() != '']
    if metadata_list[0].split(".")[-1] == "txt":
        md = [pd.read_csv(x,header=0,sep='\t') for x in metadata_list]
    if metadata_list[0].split(".")[-1] == "csv":
        md = [pd.read_csv(x,header=0,sep=',') for x in metadata_list]

if arg_dict['gene_list'] is None:
    print("No user-specified gene list available for visualization")
else:
    if isinstance(arg_dict['gene_list'], list):
        if len(arg_dict['gene_list']) > 1:
            print(str(len(arg_dict['gene_list']))+" user-specified gene lists provided.\nOnly the first user-specified gene list will be used...")
            cs_list = arg_dict['gene_list'][0].strip()
    if isinstance(arg_dict['gene_list'], str):
        cs_list = arg_dict['gene_list'].strip()
    if cs_list[-3:] == "txt":
        cs = pd.read_csv(cs_list,header=None,sep='\t')
    if cs_list[-3:] == "csv":
        cs = pd.read_csv(cs_list,header=None,sep=',')
    # format the cell-specific genes to uppercase
    cs[0] = cs[0].str.upper()

# Compile metadata columns
md = [x[[arg_dict['S'],arg_dict['T']]] for x in md]
md_final = md[0][[arg_dict['S'],arg_dict['T']]]
#for testing
md_final.iloc[33,0] = "Medial Control #2"
if len(md)>1:
    for cur_md in md[1:]:
        md_final.append(cur_md, ignore_index=True)

#Strip whitespace from sample ids
# sample_label=arg_dict['S']
md_final[arg_dict['S']] = [x.strip() for x in md_final[arg_dict['S']]]
# md_final.assign(sample_label=[x.strip() for x in md_final[arg_dict['S']]])
# df = df.assign(B=df1['E'])
# md_final[arg_dict['S']].values.strip()
#Remove duplicates samples from metadata
md_final.drop_duplicates(subset=arg_dict['S'],inplace=True)

#processing rpkms:
rpkm_final = rpkm[0]
new_columns = [x.strip() for x in rpkm_final.columns.values]; new_columns[0] = "gene_id"; rpkm_final.columns  = new_columns
rpkm_final.gene_id = rpkm_final.gene_id.str.upper()
if len(rpkm)>1:
    for cur_rpkm in rpkm[1:]:
        # cur_rpkm.head()
        new_columns = [x.strip() for x in cur_rpkm.columns.values]; new_columns[0] = "gene_id"; cur_rpkm.columns  = new_columns
        cur_rpkm.gene_id = cur_rpkm.gene_id.str.upper()
        rpkm_final = rpkm_final.merge(cur_rpkm, on="gene_id",how='inner')

print("Expression dataset dimensions:"+str(rpkm_final.shape))
#remove duplicated genes and set index
rpkm_final.drop_duplicates(subset="gene_id",keep="first",inplace=True)
rpkm_final.set_index("gene_id",inplace=True)
# print("Expression dataset dimensions:"+str(rpkm_final.shape))
#remove from rpkm duplicated columns and columns that are not in metadata
rpkm_final = rpkm_final.loc[:,~rpkm_final.columns.duplicated()]
# print("Expression dataset dimensions:"+str(rpkm_final.shape))
#remove columns that are not in metadata
rpkm_selection = rpkm_final.columns.isin(md_final[arg_dict['S']])
intersected_samples = rpkm_final.columns[rpkm_selection]
rpkm_final = rpkm_final.loc[:,rpkm_selection]
print("Final expression dataset dimensions:"+str(rpkm_final.shape))
#Check overlap of metadata sample ids with those of the expression data
# reorder and resize metadata to match rpkm
md_final.set_index(arg_dict['S'],inplace=True)
md_final.loc[intersected_samples]
md_final[arg_dict['S']] = md_final.index.values
# md_final.reset_index(inplace=True)

# Identify genes with acceptable expression values
rpkm_mean = rpkm_final.mean(axis=1)
rpkm_rmlow = rpkm_final.loc[rpkm_mean>1,:]
rpkm_sd = rpkm_rmlow.std(axis=1).sort_values(ascending=False)
# rpkm_sd
rpkm_rmlow = rpkm_rmlow.loc[rpkm_sd.index.values,:]

# Find overlaps betweeen cell-type specific genes and genes in RPKM
# gene ids in both lists are already in uppercase
if arg_dict['gene_list'] is not None:
    cs_intersect = cs[cs[0].isin(rpkm_rmlow.index.values)]
    cs_intersect.reset_index(drop=True, inplace=True)
    if cs_intersect.shape[1] > 1:
        # Assign colors
        lut = dict(zip(cs_intersect[1].unique(), sns.color_palette("hls", len(cs[1].unique()))))
        row_colors = cs_intersect[1].map(lut)
        row_colors.index = cs_intersect[0]
        # row_colors.index.name = "Cell specific\ngenes"

    rpkm_cs = rpkm_rmlow.loc[cs_intersect[0].values,:]
    print("Count of user-specified genes: "+str(cs.shape[0]))
    print("Count of user-specified genes present in expression files: "+str(cs_intersect.shape[0]))

print("Success! Now creating graphs...")

# ## heatmaps
# # do one heatmap with all genes
# # do 2nd heatmap with top most variable genes
# # do 3rd heatmap with user-specified list of cell-type specific genes
# # do all versions of heatmaps with scaled presentations and basic RPKM (or log RPKM)

#Scale features for PCA, all genes
x = rpkm_rmlow.transpose() # Separating out the target
y = md_final[arg_dict['T']] # Standardizing the features
x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
principalDf.set_index(rpkm_rmlow.columns,inplace=True)
finalDf = pd.concat([principalDf, y], axis = 1)

pca_tot = finalDf.copy()
pca_tot["id"]=finalDf.index
new_columns = ["pc1","pc2","celltype","id"]; pca_tot.columns  = new_columns

#Scale features for PCA, top 500 most variable genes
x = rpkm_rmlow.iloc[0:500,:].transpose() # Separating out the target
y = md_final[arg_dict['T']] # Standardizing the features
x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
principalDf.set_index(rpkm_rmlow.columns,inplace=True)
finalDf = pd.concat([principalDf, y], axis = 1)

pca_500 = finalDf.copy()
pca_500["id"]=finalDf.index
new_columns = ["pc1","pc2","celltype","id"]; pca_500.columns  = new_columns

heat_500 = rpkm_rmlow.iloc[0:500,:]
heat_500_log10 = np.log10(rpkm_rmlow.iloc[0:500,:]+1)
heat_tot_log10 = np.log10(rpkm_rmlow+1)
if arg_dict['gene_list'] is not None:
    heat_cs_log10 = np.log10(rpkm_cs+1)

# # improvement - plot sample name next to (or overlapping) points
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], str(point['val']),fontsize='small')

with PdfPages('multipage_pdf.pdf') as pdf:
    # scatter 1
    sns.set_style("whitegrid")
    plt.figure(figsize=(7, 7))
    ax=sns.scatterplot(data=pca_tot,x="pc1",y="pc2",hue="celltype")
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('Preliminary PCA for PCA generator, all genes', fontsize = 20)
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # scatter 2
    plt.figure(figsize=(7, 7))
    ax=sns.scatterplot(data=pca_tot,x="pc1",y="pc2",hue="celltype")
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('Preliminary PCA for PCA generator, all genes', fontsize = 20)
    label_point(pca_tot.pc1,pca_tot.pc2,pca_tot.id,plt.gca())
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # scatter 3
    plt.figure(figsize=(7, 7))
    ax=sns.scatterplot(data=pca_500,x="pc1",y="pc2",hue="celltype")
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('Preliminary PCA for PCA generator\nTop 500 most variable genes', fontsize = 20)
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # scatter 4
    plt.figure(figsize=(7, 7))
    ax=sns.scatterplot(data=pca_500,x="pc1",y="pc2",hue="celltype")
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('Preliminary PCA for PCA generator\nTop 500 most variable genes', fontsize = 20)
    label_point(pca_500.pc1,pca_500.pc2,pca_500.id,plt.gca())
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # heatmap 1
    sns.set_context(font_scale=0.5)  
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(heat_500, xticklabels=True, yticklabels=False,linewidths=0, 
        linecolor='white',cbar_kws={'label': 'RPKM'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,horizontalalignment='right',fontsize='x-small') #
    ax.set_xlabel('Samples', fontsize = 15)
    ax.set_ylabel('Genes', fontsize = 15)
    ax.set_title('Heatmap, top 500 most variable genes', fontsize = 20)
    ax.figure.subplots_adjust(bottom = 0.3)
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # heatmap 2
    sns.set_context(font_scale=0.2)  
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(heat_500_log10, xticklabels=True, yticklabels=False,linewidths=0, 
        linecolor='white',cbar_kws={'label': 'log10 RPKM'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,horizontalalignment='right',fontsize='x-small') #
    ax.set_xlabel('Samples', fontsize = 15)
    ax.set_ylabel('Genes', fontsize = 15)
    ax.set_title('Heatmap, top 500 most variable genes\nlog10 transformed RPKM', fontsize = 20)
    ax.figure.subplots_adjust(bottom = 0.3)
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # heatmap 3
    sns.set_context(font_scale=0.2)  
    plt.figure(figsize=(7, 7))
    # Assign colors
    lut = dict(zip(pca_500.celltype.unique(), sns.color_palette("hls", len(pca_500.celltype.unique()))))
    col_colors = pca_500.celltype.map(lut)
    ax = sns.clustermap(heat_500_log10, col_cluster=True, xticklabels=True, 
            yticklabels=False, cbar_kws={'label': 'log10 RPKM'},col_colors=col_colors)
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=45,horizontalalignment='right',fontsize='x-small') # For x axis
    ax.ax_heatmap.set_xlabel("Samples", fontsize = 15)
    ax.ax_heatmap.set_ylabel("Top 500 most variable genes", fontsize = 15)
    ax.ax_heatmap.annotate("Clustering of log10-transformed RPKM\nTop 500 most variable genes", xy=(0.05, 1.2), xycoords='axes fraction',
        fontsize='x-large')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    if arg_dict['gene_list'] is not None and cs_intersect.shape[0]>2:
        # heatmap 4
        sns.set_context(font_scale=0.2)  
        plt.figure(figsize=(7, 7))
        # Assign colors
        lut = dict(zip(pca_500.celltype.unique(), sns.color_palette("hls", len(pca_500.celltype.unique()))))
        col_colors = pca_500.celltype.map(lut)
        ax = sns.clustermap(heat_cs_log10, col_cluster=True,row_cluster=True, xticklabels=True, 
                yticklabels=True, cbar_kws={'label': 'log10 RPKM'},col_colors=col_colors,row_colors=row_colors)
        plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=45,horizontalalignment='right',fontsize='x-small') # For x axis
        plt.setp(ax.ax_heatmap.get_yticklabels(), rotation=0,horizontalalignment='left',fontsize='x-small') # For x axis
        ax.ax_heatmap.set_xlabel("Samples", fontsize = 15)
        ax.ax_heatmap.set_ylabel("Cell-specific genes", fontsize = 15)
        # ax.ax_heatmap.text(6.5, 0.95, "Clustering of log10-transformed RPKM", fontsize=12,
        #     verticalalignment='bottom', horizontalalignment='right')
        ax.ax_heatmap.annotate("Clustering of log10-transformed RPKM\nCell-type specific genes", xy=(0.05, 1.2), xycoords='axes fraction',
            fontsize='x-large')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

