print("This program is used for generating tSNE analysis of sequences")
print("V1.0. Code Compiled by: Varun CN, NIMHANS")
print("Require a multifasta file that contains all the sequences. Can be resource intensive depending on number of sequences")

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

print("tSNE analysis using Sklearn for WGS: V1.0, Varun CN, NIMHANS")

# Function to calculate sequence similarities
def calculate_sequence_similarities(seqs, max_seq_length):
    num_seqs = len(seqs)
    padded_seqs = []

    # Pad sequences with N (ambiguous nucleotide) to have the same length
    for seq_record in seqs:
        sequence = str(seq_record.seq)
        padding_length = max_seq_length - len(sequence)
        padded_sequence = sequence + "N" * padding_length
        padded_seqs.append(padded_sequence)

    similarities = np.zeros((num_seqs, num_seqs))

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq_i = padded_seqs[i]
            seq_j = padded_seqs[j]
            similarity = sum(1 for a, b in zip(seq_i, seq_j) if a == b) / len(seq_i)
            similarities[i, j] = similarities[j, i] = similarity

    return similarities

# Main function
def main(fasta_file_path):
    # Read the FASTA file and store sequences
    sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Find the maximum sequence length for padding
    max_seq_length = max(len(seq.seq) for seq in sequences)

    # Calculate sequence similarities with padding
    similarities = calculate_sequence_similarities(sequences, max_seq_length)

    # Perform t-SNE analysis
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(similarities)

    # Determine the optimal number of clusters using the Elbow method
    distortions = []
    max_clusters = 10  # Set a reasonable upper limit for the number of clusters
    for i in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=i, random_state=42)
        kmeans.fit(tsne_result)
        distortions.append(kmeans.inertia_)

    # Plot the Elbow curve
    #Ref: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
    plt.plot(range(1, max_clusters + 1), distortions, marker='o')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Distortion')
    plt.title('Elbow Method for Optimal Number of Clusters')
    plt.show()

    # Based on the Elbow method, enter desried number of clusters
    optimal_clusters = int(input("Optimal number of clusters as per the graph: "))

    # Perform KMeans clustering with the optimal number of clusters
    kmeans = KMeans(n_clusters=optimal_clusters, random_state=42)
    clusters = kmeans.fit_predict(tsne_result)

    # Plot the t-SNE results with different colors for each cluster
    cluster_cmap = ListedColormap(['red', 'green', 'blue'])
    plt.figure(figsize=(8, 6))
    for i in range(optimal_clusters):
        indices = np.where(clusters == i)
        plt.scatter(tsne_result[indices, 0], tsne_result[indices, 1], label=f'Cluster {i}', cmap=cluster_cmap)

    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    plt.title(f"t-SNE Analysis with {optimal_clusters} Clusters")
    plt.legend()
    plt.show()

    # Create a DataFrame to store cluster information
    cluster_data = {'Sequence ID': [seq.id for seq in sequences], 'Cluster': clusters}
    cluster_df = pd.DataFrame(cluster_data)

    # Save the cluster information to a CSV file
    cluster_df.to_csv('cluster_information.csv', index=False)

fasta_file_path = input("Provide path to multifasta file: ")

main(fasta_file_path)