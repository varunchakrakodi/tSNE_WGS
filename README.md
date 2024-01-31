# tSNE_WGS
Create a tSNE analysis for whole genome sequence data

The Python script can be used to derive a tSNE-based K-means clustering of multiple whole genome sequences for various analyses. The analysis uses the Elbow method to suggest an optimal number of clusters. The user inputs the multi-fasta file and the number of clusters based on the Elbow graph to generate a clustering graph and a .csv output file containing details of sequences from each cluster for further analysis.

**Dependencies**

numpy, pandas, biopython, sklearn, matplotlib

**Example Usage:**
python3 tSNE.py
