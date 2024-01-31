# tSNE_WGS
Create a tSNE analysis for whole genome sequence data

The Python script is based on Sklearn (https://scikit-learn.org/stable/) and can be used to derive a tSNE-based K-means clustering of multiple whole genome sequences for various analyses. The analysis uses the Elbow method (https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/) to suggest an optimal number of clusters. The user inputs the multi-fasta file and the number of clusters based on the Elbow graph to generate a clustering graph. Output includes a Elbow graph, Clustering graph and a .csv output file containing details of sequences from each cluster for further analysis.

**Dependencies**

numpy, pandas, biopython, sklearn, matplotlib

**Usage:**
python3 tSNE.py
