import os
import json
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hcluster


# Parse arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('pae', help='Path to PAE file or folder containing PAE files. If a folder is specified, script will look for files ending in *_predicted_aligned_error_v1.json.')
    parser.add_argument('-o', '--output_folder', help='Path to the output folder. Default folder is ./sublyme_domains/', default="./sublyme_domains/")
    parser.add_argument('-c', '--max_clusters', help='Maximum number of clusters defined by hierarchical clustering algorithm. (It should be a little higher than the expected number of clusters).', default="10")
    parser.add_argument('-k', '--kernel_size', help='Size of sliding window that corrects artifacts in generated clusters.', default="30")
    args = parser.parse_args()


    pae = args.pae
    output_folder = args.output_folder
    max_clusters = args.max_clusters
    kernel_size = args.kernel_size
    
    return pae, output_folder, max_clusters, kernel_size

#Loads PAE matrix from specified file.
def load_pae(filepath):
    with open(filepath) as f:
        data = json.load(f)
        
    if len(data) == 1:
        pae = np.array(data[0]["predicted_aligned_error"])
    else:
        pae = np.array(data["predicted_aligned_error"])
    return pae

# Perform hierarchical clustering to find domains
def hier_pae(pae, max_clusters=10, homogenize=True, kernel_size=30):
    pae = pae + pae.T / 2 # combine rows and columns of pae matrix to generate a symmetric matrix
    clusters = hcluster.fclusterdata(pae, max_clusters, criterion="maxclust") #hclustering
    
    if homogenize:
        clusters = homogenize_clusters(clusters, kernel_size=kernel_size) #correct artifacts and combine tiny clusters found between actual domains into linkers
        
    return clusters

# Sliding window to correct artifacts (singleton clusters within larger domains)
# as well as inter-domain linkers
def homogenize_clusters(clusters, kernel_size=10):

    half_kernel = kernel_size // 2
    for i in range(len(clusters)):
        if i >= half_kernel and i < (len(clusters) - half_kernel): #normal case
            subset = clusters[i-half_kernel:i+half_kernel]
            artifact = (not (subset == clusters[i]).all()) and (subset[0] == subset[-1])
            if artifact: #homogenize
                clusters[i-half_kernel:i+half_kernel] = subset[0]
                
        if i < half_kernel: #left border case
            subset = clusters[:i+half_kernel]
            artifact = (not (subset == clusters[i]).all()) and (subset[0] == subset[-1])
            if artifact: #homogenize
                clusters[:i+half_kernel] = subset[0]
                
        if i >= (len(clusters) - half_kernel): #right border case
            subset = clusters[i-half_kernel:]
            artifact = (not (subset == clusters[i]).all()) and (subset[0] == subset[-1])
            if artifact: #homogenize
                clusters[i-half_kernel:] = subset[0]

    #Replace clusters <kernel_size amino acids long by a single cluster number
    outliers = np.unique(clusters, return_counts=True)[0][np.unique(clusters, return_counts=True)[1] < kernel_size]
    clusters = np.where(np.isin(clusters, outliers), 42, clusters)

    #Treat edge cases. Place them in nearest cluster
    if not (clusters[:10] == clusters[0]).all():
        clusters[:10] = clusters[11]
    if not (clusters[-10:] == clusters[-1]).all():
        clusters[-10:] = clusters[-11]
                
    return clusters

#Find the positions of residues corresponding to domain boundaries
def find_delineations(clusters):
    j=clusters[0]
    delineations=[1]
    for i, clust in enumerate(clusters):
        if clust != j:
            delineations.append(i+1)
            j=clust
    delineations.append(i+1)
    return delineations


#plot pae with clusters
def plot_pae(pae, clusters, output_folder, sample_name):
    cols = map_cols(clusters, len(np.unique(clusters)))
    g = sns.clustermap(data=pae, col_colors=cols, row_cluster=False, col_cluster=False, dendrogram_ratio=0.05, cbar_pos=None, figsize=(5, 5))
    plt.savefig(os.path.join(output_folder, f"{sample_name}_pae.png"))

#map colors to clusters for plotting
def map_cols(clusters, n):
    colors = distinctipy.get_colors(n)
    
    c_dict={}
    for i, clust in enumerate(np.unique(clusters)):
        c_dict[clust] = colors[i]

    cols=[]
    for i in clusters:
        cols.append(c_dict[i])
    return cols


def sublyme(pae, output_folder, max_clusters, kernel_size):
    pattern = "_predicted_aligned_error_v1.json"

    #Check that specified path exists
    if not os.path.exists(pae_path):
        raise FileNotFoundError(f"{filepath} was not found.")
        
    #If input is a file
    if os.path.isfile(pae_path):
        sample_name = os.path.basename(pae_path).rstrip(pattern) #find sample name
        pae = load_pae_(pae_path) #load pae matrix

        clusters = hier_pae(pae, max_clusters=max_clusters, homogenize=True, kernel_size=kernel_size) #find clusters

        delineations = find_delineations(clusters) #find domain boundaries

        output = pd.DataFrame(data=[sample_name, delineations])
        output.to_csv(os.path.join(output_folder, f"{sample_name}_delineations.csv"))
        plot_pae(pae, clusters, output_folder, sample_name)
        
        

    #If input is a folder
    else:
        pae_filenames = [f for f in os.listdir(pae_path) if pattern in f] #find filenames
        all_delineations = pd.DataFrame()
        
        for file in pae_filenames:
            filepath = os.path.join(pae_path, file)
            
            sample_name = os.path.basename(filepath).rstrip(pattern) #find sample name
            pae = load_pae_(os.path.join(filepath, file)) #load pae matrix

            clusters = hier_pae(pae, thresh=10, cluster_on="both", homogenize=True, kernel_size=30) #find clusters

            delineations = find_delineations(clusters) #find domain boundaries

            plot_pae(pae, clusters, output_folder, sample_name)
            
            all_delineations = pd.concat(all_delineations, pd.DataFrame(data=[sample_name, delineations]))
            
        all_delineations.to_csv(os.path.join(output_folder, "delineations.csv"))
            


if __name__=='__main__':
    pae_path, output_folder, max_clusters, kernel_size = argparse()
    
    sublyme(pae_path, output_folder, max_clusters, kernel_size)