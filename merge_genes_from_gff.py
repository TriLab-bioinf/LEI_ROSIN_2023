#!/Users/lorenziha/opt/anaconda3/bin/python

import pandas as pd
import re
import argparse

# Process parameters

# Initialize parser
parser = argparse.ArgumentParser(
    description="This proram cluster genes located on the same strand of a gff file that overlap more than a defined X% (default = 0.5 => 50%) of the smallest gene\n\n"
)

# Adding argument
parser.add_argument(
    "-g",
    "--gff",
    help="gff file",
    type=str,
    required=True,
)

parser.add_argument(
    "-o",
    "--overlap",
    help="Percent overlap of smallest gene [0 to 1, default = 0.5]",
    type=float,
    required=False,
    default=0.5,
)

args = parser.parse_args()
input_file = args.gff
PERCENT_OVERLAP = args.overlap
CLUSTER_ID = 0

# Functions
def merge_genes(df):
    allLines = df.to_string(header=False, index=False, index_names=False).split('\n')
    #allLines = ['\t'.join(ele.split(' ')) for ele in allLines]
    #df.sort_values(ascending=True, kind='mergesort', by=[3])
    arr5 = list(df[3])
    arr3 = list(df[4])
    index = 0
    # print(allLines)
    # print(f"index={index}",arr5[0],arr3[0])
    arr = [[arr5[0],arr3[0]]]
    genes = [[re.search(r'ID="(\S*)";', str(allLines[0])).group(1)]]
    lines = [[allLines[0]]]
    for i in range(1, len(arr5)):
        # If this is not first Interval and overlaps
        # with the previous one, Merge previous and
        # current Intervals
        gene = re.search(r'ID="(\S*)";', str(allLines[i]))
        #print(f'gene={gene.group(1)}')
        if (arr[index][1] >= arr5[i]):
            overlap_len = arr[index][1] - arr5[i]
            upstream_overlap = (overlap_len / (arr[index][1] - arr[index][0]))
            downstream_overlap = (overlap_len / (arr3[i] - arr5[i]))
            #print(f'Upstream = {overlap_len}', f'Downstream = {downstream_overlap}')
            if (upstream_overlap > PERCENT_OVERLAP or downstream_overlap > PERCENT_OVERLAP):
                arr[index][1] = max(arr[index][1], arr3[i])
                #print(f'index={index}')
                #print(genes)
                #print(f'genes={len(genes)}, arr={len(arr)}')
                genes[-1].append(gene.group(1))
                lines[-1].append(allLines[i])
            else:
                index = index + 1
                arr.append([arr5[i],arr3[i]])
                genes.append([gene.group(1)])
                lines.append([allLines[i]])
        else:
            index = index + 1
            arr.append([arr5[i],arr3[i]])
            genes.append([gene.group(1)])
            lines.append([allLines[i]])
    #print("The Merged Intervals are :", end=" ")
    global CLUSTER_ID
    for i in range(len(genes)):
        CLUSTER_ID = CLUSTER_ID + 1
        for j in range(len(genes[i])):
            #print(genes[i][j], CLUSTER_ID)
            my_line = lines[i][j].split(maxsplit=8)
            my_line.append(str(CLUSTER_ID))
            print('\t'.join(my_line))
            #print("\t".join(lines[i][j].split(maxsplit=7)), CLUSTER_ID)
            #print(lines[i][j], CLUSTER_ID)


# Import gff file
gff = pd.read_table(filepath_or_buffer=input_file, header=None)
genes_gff = gff[gff[2] == "gene"]
genes_gff = genes_gff.replace('; ',';')


# Iterate through each contig and strands
for i in genes_gff[0].unique():
    # print(i)
    my_df_fwd = genes_gff[(genes_gff[0] == i) & (genes_gff[6] == '+')]
    if (len(my_df_fwd) > 0):
        # print(my_df_fwd)
        merge_genes(my_df_fwd)
    my_df_rev = genes_gff[(genes_gff[0] == i) & (genes_gff[6] == '-')]
    if (len(my_df_rev) > 0):
        merge_genes(my_df_rev)