import os
import json
import argparse
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as hcluster


# Parse arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('pae_path', help='Path to PAE file or folder containing PAE files. If a folder is specified, script will look for files ending in *.json.')
    parser.add_argument('-o', '--output_file', help='Path to the output file. Default file is ./spaed_predictions.csv', default="./spaed_predictions.csv")
    parser.add_argument('--MAX_CLUSTERS', type=str, help='Maximum number of clusters initially generated by hierarchical clustering (default "dynamic"). By default will use a value equal to len(protein) // 10. Can also accept any integer.', default="dynamic")
    parser.add_argument('--MIN_DOMAIN_SIZE', type=int, help='Minimum size a domain can have (default 30 residues).', default=30)
    parser.add_argument('--MIN_DISORDERED_SIZE', type=int, help='Minimum size a disordered region can be to be interesting enough to separate from the domain it is next to (default 20 residues).', default=20)
    parser.add_argument('--PAE_SCORE_CUTOFF', type=float, help='Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions. Residues with PAE score < PAE_SCORE_CUTOFF are considered close together. (default 5)', default=5)
    parser.add_argument('--FREQ_DISORDERED', type=int, help='For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values <MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted). (default 6)', default=6)
    parser.add_argument('--PROP_DISORDERED', type=float, help="Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict the region as disordered.", default=0.80)
    parser.add_argument('--FREQ_LINKER', type=int, help='For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered as part of the linker. Values < MIN_DOMAIN_SIZE are logical as they are less than the expected size of the nearest domain. Increasing leads to a more leniant assignment of residues as part of the linker. (default 15)', default=15)
    args = parser.parse_args()


    pae_path = args.pae_path
    output_file = args.output_file
    MAX_CLUSTERS = args.MAX_CLUSTERS
    if MAX_CLUSTERS != "dynamic": MAX_CLUSTERS = int(MAX_CLUSTERS)
    MIN_DOMAIN_SIZE = args.MIN_DOMAIN_SIZE
    MIN_DISORDERED_SIZE = args.MIN_DISORDERED_SIZE
    PAE_SCORE_CUTOFF = args.PAE_SCORE_CUTOFF
    FREQ_DISORDERED = args.FREQ_DISORDERED
    PROP_DISORDERED = args.PROP_DISORDERED
    FREQ_LINKER = args.FREQ_LINKER

    return pae_path, output_file, MAX_CLUSTERS, MIN_DOMAIN_SIZE, MIN_DISORDERED_SIZE, PAE_SCORE_CUTOFF, FREQ_DISORDERED, PROP_DISORDERED, FREQ_LINKER

"""
User parameters:
- MAX_CLUSTERS: Maximum number of clusters initially generated by hierarchical clustering
- MIN_DOMAIN_SIZE: Minimum size a domain can have
- MIN_DISORDERED_SIZE: Minimum size a disordered region can be to be interesting enough to separate from the domain it is next to
- PAE_SCORE_CUTOFF: Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions
- FREQ_DISORDERED: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values <MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted)
- PROP_DISORDERED: Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality a little more conservative); more than 50% of residues must be associated to a disordered region to consider this region as disordered.
- FREQ_LINKER: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered as part of the linker. Values < MIN_DOMAIN_SIZE are logical as they are less than the expected size of the nearest domain. Increasing leads to a more leniant assignment of residues as part of the linker.

Development parameters (to make debugging easier):
- form: whether to perform first step after clustering (predicting which clusters are part of domains)
- ends: whether to perform end correction step
- artifacts: whether to perform artifact correction step
- linkers: whether to perform linker adjustment step
- plot: whether to plot pae matrix with predicted domains

Returns:
- which domain each residue corresponds to (list of ints)
- Number of predicted domains (excluding linkers and disordered regions)
"""
def spaed_(pae, MAX_CLUSTERS=40, MIN_DOMAIN_SIZE=30, MIN_DISORDERED_SIZE=20, PAE_SCORE_CUTOFF=5, FREQ_DISORDERED=6,
             PROP_DISORDERED=0.80, FREQ_LINKER=15, form=True, ends=True, artifacts=True, linkers=True, plot=False):

    if MAX_CLUSTERS=="dynamic":
        clusters = hcluster.fclusterdata(pae, len(pae)//10, criterion="maxclust") #perform clustering
    else:
        clusters = hcluster.fclusterdata(pae, MAX_CLUSTERS, criterion="maxclust") #perform clustering

    if form: #identify domains
        clusters = form_domains(pae, clusters, MIN_DOMAIN_SIZE, PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE, FREQ_DISORDERED, PROP_DISORDERED)

    if artifacts: #correct artifacts
        clusters = correct_artifacts(clusters)

    if ends: #correct the ends of sequences and detect disordered regions
        clusters = correct_ends(pae, clusters, PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE, FREQ_DISORDERED, PROP_DISORDERED, MIN_DOMAIN_SIZE)

    if linkers: #adjust linkers
        adjust_linkers(pae, clusters, PAE_SCORE_CUTOFF, MIN_DOMAIN_SIZE, FREQ_DISORDERED, FREQ_LINKER)

    if plot: #plot pae matrix with predicted domains
        cols = map_cols(clusters, len(np.unique(clusters)), form=form)
        g = sns.clustermap(data=pae, col_colors=cols, row_cluster=False, col_cluster=False, dendrogram_ratio=0.05, cbar_pos=None, figsize=(5, 5))
        plt.show()

    clusters = norm_cluster_numbers(clusters)

    num_domains = len(set(clusters.unique()) - set([-1, -2]))

    return clusters, num_domains


"""
Function to load PAE matrix.

- filepath: path to PAE json file

Returns: normalized PAE matrix
"""
def load_pae(filepath):
    with open(filepath) as f:
        data = json.load(f)

    if len(data) == 1:
        pae = np.array(data[0]["predicted_aligned_error"])
    else:
        pae = np.array(data["predicted_aligned_error"])

    pae = pae + pae.T / 2 #make matrix symmetrical

    return pae


"""
Function to predict which clusters are domains.

- pae: normalized pae matrix
- clusters: clusters each residue were assigned to at this step
- MIN_DOMAIN_SIZE: Minimum size a domain can have
- PAE_SCORE_CUTOFF: Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions
- MIN_DISORDERED_SIZE: Minimum size a disordered region can be to be interesting enough to separate from the domain it is next to
- FREQ_DISORDERED: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values < MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted).
- PROP_DISORDERED: Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality set a little more conservative)1~value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality set a little more conservative); more than 50% of residues must be associated to a disordered region to consider this region as disordered.

Returns: adjusted clusters each residue are assigned to after treatment
"""
def form_domains(pae, clusters, MIN_DOMAIN_SIZE, PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE, FREQ_DISORDERED, PROP_DISORDERED):
    clusters = pd.Series(clusters)
    counts = clusters.value_counts()

    domain_ids = counts.loc[counts >= MIN_DOMAIN_SIZE-5].index #threshold for size of domains

    clusters.loc[~clusters.isin(domain_ids)] = -1 #Distinguish domains from linkers, extremities and artifacts

    #if two clusters are predicted as domains, but are overlapping (ex. aaaabbbbaaaabbbbccccccc) they are merged as 1 domain.
    indices = []
    for i, c_1 in enumerate(domain_ids):
        for j, c_2 in enumerate(domain_ids):
            if ((i != j) & (i < j)):
                i_min = clusters[clusters == c_1].index[0]; i_max = clusters[clusters == c_1].index[-1]
                j_min = clusters[clusters == c_2].index[0]; j_max = clusters[clusters == c_2].index[-1]

                if not (((i_min < j_min) & (i_max < j_min)) | ((j_min < i_min) & (j_max < i_min))):
                    clusters[clusters == c_1] = c_2
                    domain_ids = domain_ids.drop(c_1)
    for c in domain_ids:
        c_domain = clusters.loc[clusters == c].index
        disordered, _ = detect_disordered(pae, (c_domain[0], c_domain[-1]), clusters, "domain", PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE,
                                                 FREQ_DISORDERED, PROP_DISORDERED) #detect domains that are in fact disordered
        if disordered:
            clusters[c_domain[0] : c_domain[-1] + 1] = -1

    return clusters


"""
Function to correct the end clusters (either as part of previous domain or as a disordered region).

- pae: normalized pae matrix
- clusters: clusters each residue were assigned to at this step
- PAE_SCORE_CUTOFF: Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions
- MIN_DISORDERED_SIZE: Minimum size a disordered region can be to be interesting enough to separate from the domain it is next to
- FREQ_DISORDERED: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values < MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted).
- PROP_DISORDERED: Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality set a little more conservative)1~value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality set a little more conservative); more than 50% of residues must be associated to a disordered region to consider this region as disordered.
- MIN_DOMAIN_SIZE: Minimum size a domain can have

Returns: adjusted clusters each residue are assigned to after treatment
"""
def correct_ends(pae, clusters, PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE, FREQ_DISORDERED, PROP_DISORDERED, MIN_DOMAIN_SIZE):
    start_first = clusters.loc[clusters != -1].index.min()
    end_last = clusters.loc[clusters != -1].index.max()

    disordered = False
    #Beginnning of sequence
    if start_first > MIN_DISORDERED_SIZE: #Around 20-25
        disordered, clusters = detect_disordered(pae, (0, start_first), clusters, "start", PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE,
                                                 FREQ_DISORDERED, PROP_DISORDERED) #detect disordered regions
    if (not disordered): clusters[:start_first] = clusters[start_first]

    disordered = False
    #End of sequence
    if len(clusters)-end_last > MIN_DISORDERED_SIZE: #Around 20-25
        disordered, clusters = detect_disordered(pae, (end_last+1, len(clusters)), clusters, "end", PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE,
                                                  FREQ_DISORDERED, PROP_DISORDERED) #detect disordered regions
    if (not disordered): clusters[end_last+1:] = clusters[end_last]

    return clusters


"""
Function to check for disordered regions within pae matrix (at specified location: beginning or end of sequence)

- pae: normalized pae matrix
- subset_indices: position of residues in region that is being evaluated
- clusters: clusters each residue were assigned to at this step
- region: "start" or "end"
- PAE_SCORE_CUTOFF: Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions
- MIN_DISORDERED_SIZE: Minimum size a disordered region can be to be interesting enough to separate from the domain it is next to
- FREQ_DISORDERED: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values < MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted).
- PROP_DISORDERED: Proportion of residues in a given region that must meet FREQ_DISORDERED criteria to be considered a disordered region. The greater the value, the stricter the criteria to predict a region as disordered. By definition must be higher than 50% (in reality set a little more conservative); more than 50% of residues must be associated to a disordered region to consider this region as disordered.

Returns:
- True if region is disordered,
- adjusted clusters each residue are assigned to after treatment
"""
def detect_disordered(pae, subset_indices, clusters, region, PAE_SCORE_CUTOFF, MIN_DISORDERED_SIZE, FREQ_DISORDERED, PROP_DISORDERED):
    pae_counts = pd.Series((pae < PAE_SCORE_CUTOFF).sum(axis=1)).loc[subset_indices[0]:subset_indices[1]+1]

    if (pae_counts < FREQ_DISORDERED).sum() == 0: return False, clusters #no disordered residues in region

    #Find proportion of region that is disordered (PAE<6)
    if region == "start":
        l_border=0; r_border=pae_counts[pae_counts < FREQ_DISORDERED].index[-1]
        prop = (pae_counts.loc[l_border:r_border+1] <= FREQ_DISORDERED).sum() / (r_border+1 - l_border)

    elif region == "end":
        l_border=pae_counts[pae_counts < FREQ_DISORDERED].index[0]; r_border=len(clusters)-1
        prop = (pae_counts.loc[l_border:r_border+1] <= FREQ_DISORDERED).sum() / (r_border+1 - l_border)

    elif region == "domain":
        l_border=pae_counts.index[0]; r_border=pae_counts.index[-1]
        prop = (pae_counts <= FREQ_DISORDERED).sum() / len(pae_counts)

    if (r_border+1 - l_border) < MIN_DISORDERED_SIZE: return False, clusters #too short to care about it being disordered

    disordered = prop > PROP_DISORDERED #more than half of residues in this region are disordered
    if disordered & ((region=="start") | (region=="end")):
        clusters[l_border:r_border+1] = -2 #ID for disordered regions

    return disordered, clusters


"""
Function to correct artifacts (little errors in the middle of predicted domains).

- clusters: clusters each residue were assigned to at this step

Returns: clusters each residue are assigned to after this step
"""
def correct_artifacts(clusters):
    correction_value = None
    inds_to_correct = []

    for i, val in enumerate(clusters):
        if not val in [-1, -2, correction_value]:
            correction_value = val
            inds_to_correct = []
        if (val == -1) and (correction_value != None):
            inds_to_correct += [i]
        if (val == correction_value) & (inds_to_correct!=[]):
            clusters.loc[inds_to_correct] = correction_value
    return clusters


"""
Function to find the indices of all linkers (and to split the different linkers)

- clusters: clusters each residue were assigned to at this step

Returns: the indices of each linker (list of lists)
"""
def get_linker_indices(clusters):
    linkers = clusters.loc[clusters == -1].index
    if len(linkers) == 0: return []

    group_indices = []
    indices=[linkers[0]]
    for i in range(1, len(linkers)):
        if linkers[i] == linkers[i - 1] + 1:
            indices.append(linkers[i])
        else:
            group_indices.append(indices)
            indices = [linkers[i]]
    group_indices.append(indices)

    return group_indices

"""
Function to adjust boundaries of 1 linker.

- pae_counts: counts of rows in pae matrix that have a score <= PAE_SCORE_CUTOFF
- clusters: clusters each residue were assigned to at this step
- linker: indices (position in sequence) of residues that are part of a linker
- FREQ_LINKER: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered as part of the linker. Values < MIN_DOMAIN_SIZE are logical as they are less than the expected size of the nearest domain. Increasing leads to a more leniant assignment of residues as part of the linker.

Returns: adjusted clusters each residue are assigned to after treatment
"""
def adjust_linker_(pae_counts, clusters, linker, FREQ_LINKER):
    l_border = linker[0]; r_border=linker[-1]
    thresh  = pae_counts[l_border : r_border+1].min() + FREQ_LINKER #Thresh depends on linker and level of packing in PAE matrix

    #if linker is next to a disorded region, it cannot be a linker (either add to linker or to nearest domain)
    if (clusters[l_border-1] == -2): #left of linker
        for i in linker:
            if pae_counts[i] < FREQ_LINKER: clusters[i] = -2
            else: clusters[i:r_border+1] = clusters[r_border+1]; return clusters
        return clusters

    if (clusters[r_border+1] == -2): #right of linker
        for i in reversed(linker):
            if pae_counts[i] < FREQ_LINKER: clusters[i] = -2
            else: clusters[l_border:i+1] = clusters[l_border-1]; return clusters
        return clusters

    linker_tr = pae_counts[l_border: r_border+1] <= thresh #find "core" of linker
    l_border_new = pd.Series(linker_tr)[linker_tr].index[0]; r_border_new = pd.Series(linker_tr)[linker_tr].index[-1]

    # Left boundary of linker
    if (l_border_new == 0) & (pae_counts[l_border] <= thresh):
        extend = True #extend border
        i=0
        while (extend == True) & (l_border-i != 0):
            clusters[l_border-i] = -1
            i+=1
            extend = pae_counts[l_border - i] <= thresh
    else:
        clusters[l_border:l_border+l_border_new] = clusters[l_border-1] #chop

    #right boundary of linker
    if r_border_new == len(linker_tr)-1:
        extend = True #extend border
        i=0
        while (extend == True) & (r_border+i < len(clusters)-1):
            clusters[r_border+i] = -1
            i += 1
            extend = pae_counts[r_border+i] <= thresh
            if (r_border+i == len(clusters)-1): clusters[r_border+i] = -1
    else:
        clusters[l_border+r_border_new+1:r_border+1] = clusters[r_border+1] #chop

    return clusters


"""
Function to adjust the boundaries of linkers

- pae: normalized pae matrix
- clusters: clusters each residue were assigned to at this step
- PAE_SCORE_CUTOFF: Cutoff on the PAE score used to make adjustments to predicted domains/linkers/disordered regions
- MIN_DOMAIN_SIZE: Minimum size a domain can have
- FREQ_DISORDERED: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered "not part of a domain". Values < MIN_DOMAIN_SIZE are logical, but as it increases, the more leniant the algorithm becomes to disordered regions (more will be predicted).
- FREQ_LINKER: For a given residue in the PAE matrix, frequency of residues that can align to it with a low PAE score and still be considered as part of the linker. Values < MIN_DOMAIN_SIZE are logical as they are less than the expected size of the nearest domain. Increasing leads to a more leniant assignment of residues as part of the linker.

- Returns: adjusted clusters each residue are assigned to after treatment
"""
def adjust_linkers(pae, clusters, PAE_SCORE_CUTOFF, MIN_DOMAIN_SIZE, FREQ_DISORDERED, FREQ_LINKER):
    if (clusters==-1).sum() == 0: return clusters
    linkers = get_linker_indices(clusters)

    pae_counts = (pae < PAE_SCORE_CUTOFF).sum(axis=1)

    for linker in linkers:
        clusters = adjust_linker_(pae_counts, clusters, linker, FREQ_LINKER)

    #check if long linkers are actually domains that got missed
    linkers = get_linker_indices(clusters)
    for linker in linkers:
        if len(linker) > MIN_DOMAIN_SIZE-10: #probably a domain
            prop = (pae_counts[linker] < FREQ_DISORDERED).sum() / len(linker) #proportion of disordered residues
            domain = prop < 0.50 #most probably a domain surrounded by two linkers; at least half of residues are ordered

            if domain:
                #If linker is at end of sequence because of previous correction, set it as a disordered region
                if (pae_counts[linker] <= FREQ_LINKER).all(): clusters[linker] = -2; return clusters

                #find actual linkers on either side of domain
                dom_start = -1; dom_end = -1
                for i in linker:
                    if pae_counts[i] > FREQ_LINKER:
                        dom_start = i; break
                for i in reversed(linker):
                    if pae_counts[i] > FREQ_LINKER:
                        dom_end = i; break
                clusters[dom_start : dom_end+1] = 102 #assign a new cluster number to the middle region

    return clusters


"""
Function to normalize (domains are assigned 0, 1, 2...) the final cluster numbers that are assigned.

-clusters: clusters each residue were assigned to at this step

Returns: the normalized cluster numbers assigned to each residue
"""
def norm_cluster_numbers(clusters):
    i=1000

    uni_indexes = np.unique(clusters, return_index=True)[1] #np.unique without reordering
    unique_clusters = [clusters[index] for index in sorted(uni_indexes)]
    for value in unique_clusters:
        if value == -1:
            pass
        elif value == -2:
            pass
        else:
            clusters[clusters == value] = i
            i+=1
    clusters=clusters % 1000
    clusters[clusters == 998] = -2
    clusters[clusters == 999] = -1
    return clusters


#Function to map colors to clusters to plot pae matrix with clusters
def map_cols(clusters, n, form=True):
    if not form:
        colors = distinctipy.get_colors(n+2)
    else:
        colors = [[0.8392156862745098, 0.15294117647058825, 0.1568627450980392], [0.4980392156862745, 0.4980392156862745, 0.4980392156862745],
                [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], [1.0, 0.4980392156862745, 0.054901960784313725],
                [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], [0.5803921568627451, 0.403921568627451, 0.7411764705882353],
                [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], [0.8901960784313725, 0.4666666666666667, 0.7607843137254902],
                [0.7372549019607844, 0.7411764705882353, 0.13333333333333333], [0.09019607843137255, 0.7450980392156863, 0.8117647058823529]]

    c_dict={}; c=2
    uni_indexes = np.unique(clusters, return_index=True)[1]
    unique_clusters = [clusters[index] for index in sorted(uni_indexes)]
    for clust in unique_clusters:
        if clust == -1:
            c_dict[clust] = colors[1]
        elif clust == -2:
            c_dict[clust] = colors[0]
        else:
            c_dict[clust] = colors[c]
            c+=1

    cols=[]
    for i in clusters:
        cols.append(c_dict[i])
    return cols


def get_delineations(clusters):
    delin = pd.DataFrame(["", "", ""], index=["domains", "linkers", "disordered"]).T

    clust = clusters[0]
    start = 0; end = 0
    for i in range(len(clusters)):
        if clust != clusters[i]:
            end = i-1

            if clust == -1: delin.loc[0, "linkers"] = delin.loc[0, "linkers"] + f";{start+1}-{end+1}"
            elif clust == -2: delin.loc[0, "disordered"] = delin.loc[0, "disordered"] + f";{start+1}-{end+1}"
            else: delin.loc[0, "domains"] = delin.loc[0, "domains"] + f";{start+1}-{end+1}"

            start = i
            clust = clusters[i]
    end = len(clusters) - 1
    if clust == -1: delin.loc[0, "linkers"] = delin.loc[0, "linkers"] + f";{start+1}-{end+1}"
    elif clust == -2: delin.loc[0, "disordered"] = delin.loc[0, "disordered"] + f";{start+1}-{end+1}"
    else: delin.loc[0, "domains"] = delin.loc[0, "domains"] + f";{start+1}-{end+1}"

    delin.linkers = delin.linkers.str.lstrip(";")
    delin.disordered = delin.disordered.str.lstrip(";")
    delin.domains = delin.domains.str.lstrip(";")

    return delin


"""
Function that launches the whole pipeline.
"""
def spaed(pae_path, output_file="./spaed_predictions.csv", MAX_CLUSTERS="dynamic", MIN_DOMAIN_SIZE=30, MIN_DISORDERED_SIZE=20, PAE_SCORE_CUTOFF=5, FREQ_DISORDERED=6, PROP_DISORDERED=0.80, FREQ_LINKER=15):
    pattern = "_predicted_aligned_error_v1"
    all_delineations = pd.DataFrame(columns=["length", "# domains", "domains", "linkers", "disordered"]).astype(object)
    print("Start")
    #Check that specified path to input exists
    if not os.path.exists(pae_path):
        raise FileNotFoundError(f"{pae} was not found.")

    #Check that specified path to output folder exists
    if not os.path.exists(os.path.dirname(f"./{output_file}")):
        raise FileNotFoundError(f"{output_file} was not found.")

    if os.path.isfile(pae_path):
        sample_name = os.path.basename(pae_path).replace(".json", "").replace(pattern, "")

        try:
            pae = load_pae(pae_path)
        except:
            print(f"File not formatted correctly: {pae_path}")

        clusters, num_domains = spaed_(pae, MAX_CLUSTERS, MIN_DOMAIN_SIZE, MIN_DISORDERED_SIZE, PAE_SCORE_CUTOFF, FREQ_DISORDERED, PROP_DISORDERED, FREQ_LINKER)
        delin = get_delineations(clusters)
        all_delineations.loc[sample_name, ["domains", "linkers", "disordered"]] = delin.iloc[0]
        all_delineations.loc[sample_name, "length"] = len(pae)
        all_delineations.loc[sample_name, "# domains"] = num_domains

    else:
        pae_filenames = [f for f in os.listdir(pae_path) if f.endswith(".json")] #find filenames

        for file in pae_filenames:
            filepath = os.path.join(pae_path, file)
            sample_name = os.path.basename(filepath).replace(".json", "").replace(pattern, "") #find sample name
            print(sample_name)
            try:
                pae = load_pae(filepath)
            except:
                print(f"File not properly formatted: {filepath}")
                continue

            try:
                clusters, num_domains = spaed_(pae, MAX_CLUSTERS, MIN_DOMAIN_SIZE, MIN_DISORDERED_SIZE, PAE_SCORE_CUTOFF, FREQ_DISORDERED, PROP_DISORDERED, FREQ_LINKER)
            except:
                print(f"Error with {sample_name}. File was skipped, please investigate.")
                continue
            delin = get_delineations(clusters)
            all_delineations.loc[sample_name, ["domains", "linkers", "disordered"]] = delin.loc[0]
            all_delineations.loc[sample_name, "length"] = len(pae)
            all_delineations.loc[sample_name, "# domains"] = num_domains

    all_delineations.to_csv(output_file)



if __name__=='__main__':
    pae_path, output_file, MAX_CLUSTERS, MIN_DOMAIN_SIZE, MIN_DISORDERED_SIZE, PAE_SCORE_CUTOFF, FREQ_DISORDERED, PROP_DISORDERED, FREQ_LINKER = parse_args()

    spaed(pae_path, output_file, MAX_CLUSTERS, MIN_DOMAIN_SIZE, MIN_DISORDERED_SIZE, PAE_SCORE_CUTOFF, FREQ_DISORDERED, PROP_DISORDERED, FREQ_LINKER)
