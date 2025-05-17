import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
from scipy.stats import rankdata
import h5py
import re


# making a dictionary to map genes to regions. Each gene will map to an
# arrray of regions
def map_genes_to_regions(all_regions_gene_ids):
    gene2region = {}
    for genes, region in zip(all_regions_gene_ids,
                             range(len(all_regions_gene_ids))):
        if len(genes) > 1:
            for gene in genes:
                if gene not in gene2region:
                    gene2region[gene] = np.array([region])
                else:
                    gene2region[gene] = np.append(gene2region[gene],
                                                  [region])
    return gene2region


# making dictionaries to map a gene to similar genes and corresponding scores
def map_genes_to_genes(metric_top1000):
    gene2genes = {}
    gene2scores = {}
    for row in metric_top1000:
        gene0 = row[0]
        gene1 = row[1]
        score = int(row[2])
        if gene0 not in gene2genes:
            gene2genes[gene0] = np.array([gene1])
            gene2scores[gene0] = np.array([score])
        else:
            gene2genes[gene0] = np.append(gene2genes[gene0], [gene1])
            gene2scores[gene0] = np.append(gene2scores[gene0], [score])
    return gene2genes, gene2scores


# converting string scores to a num_regions x num_saved_recommended_regions
# matrix
def create_sim_matrix(all_regions_gene_ids, metric_top1000, gene2region,
                      gene2genes, gene2scores, starting_region):
    genes = np.unique(all_regions_gene_ids[starting_region])  # unique
    # if there are no genes
    if len(genes) == 0:
        region_scores_to_m = np.zeros(1000, dtype=int)
        similar_regions_to_m = np.arange(1000)
        origin_gene_tracker_to_m = np.array(['nan']*1000)
        destination_gene_tracker_to_m = np.array(['nan']*1000)
    # if there is at least 1 gene
    else:
        region_scores_init = np.array([], dtype=int)
        similar_genes_init = np.array([])
        origin_gene_tracker_init = np.array([])
        for gene in genes:
            # finding similar genes/scores based on origin gene
            region_scores_init = np.append(region_scores_init,
                                           gene2scores[gene])
            similar_genes_init = np.append(similar_genes_init,
                                           gene2genes[gene])
            scores_len = len(gene2scores[gene])
            origin_gene_tracker_init = np.append(origin_gene_tracker_init,
                                                 np.array([gene] *
                                                          scores_len))

        region_scores = np.array([], dtype=int)
        similar_regions = np.array([], dtype=int)
        origin_gene_tracker = np.array([])
        destination_gene_tracker = np.array([])
        for gene, score, origin_gene in zip(similar_genes_init,
                                            region_scores_init,
                                            origin_gene_tracker_init):
            # if gene in some region, add scores/regions to final arrays
            if gene in gene2region:
                gene_regions = np.unique(gene2region[gene])  # unique
                similar_regions = np.append(similar_regions, gene_regions)
                region_scores = np.append(region_scores,
                                          np.array([score] *
                                                   len(gene_regions)))
                origin_gene_tracker = np.append(origin_gene_tracker,
                                                np.array([origin_gene] *
                                                         len(gene_regions)))
                gr_len = len(gene_regions)
                destination_gene_tracker = np.append(destination_gene_tracker,
                                                     np.array([gene] *
                                                              gr_len))

        # sort scores
        sort_order = np.argsort(-region_scores)
        region_scores = region_scores[sort_order]
        similar_regions = similar_regions[sort_order]
        origin_gene_tracker = origin_gene_tracker[sort_order]
        destination_gene_tracker = destination_gene_tracker[sort_order]

        # finding first occurrence
        _, first_occurrence_ind = np.unique(similar_regions,
                                            return_index=True)
        first_occurrence_ind_sorted = np.sort(first_occurrence_ind)

        # values with first occurrence
        similar_regions_to_m = similar_regions[first_occurrence_ind_sorted]
        region_scores_to_m = region_scores[first_occurrence_ind_sorted]
        origin_gene_tracker_to_m = origin_gene_tracker[
                                       first_occurrence_ind_sorted
                                   ]
        destination_gene_tracker_to_m = destination_gene_tracker[
                                            first_occurrence_ind_sorted
                                        ]

        # removing starting region values, take top 1000 regions/scores
        starting_region_mask = similar_regions_to_m != starting_region
        similar_regions_to_m = similar_regions_to_m[
                                   starting_region_mask
                               ][:1000]
        region_scores_to_m = region_scores_to_m[
                                 starting_region_mask
                             ][:1000]
        origin_gene_tracker_to_m = origin_gene_tracker_to_m[
                                       starting_region_mask
                                   ][:1000]
        destination_gene_tracker_to_m = destination_gene_tracker_to_m[
                                            starting_region_mask
                                        ][:1000]

        arr_len = len(region_scores_to_m)
        if arr_len < 1000:
            region_scores_to_m = np.append(region_scores_to_m,
                                           np.zeros(1000 - arr_len,
                                                    dtype=int))
            similar_regions_to_m = np.append(similar_regions_to_m,
                                             np.array([starting_region] *
                                                      (1000-arr_len)))
            origin_gene_tracker_to_m = np.append(origin_gene_tracker_to_m,
                                                 np.array(['nan'] *
                                                          (1000 - arr_len)))
            destination_gene_tracker_to_m = np.append(
                                                destination_gene_tracker_to_m,
                                                np.array(['nan'] *
                                                         (1000 - arr_len))
                                            )

    return np.append(np.append(np.append(region_scores_to_m.astype(str),
                                         similar_regions_to_m.astype(str)),
                               origin_gene_tracker_to_m),
                     destination_gene_tracker_to_m)


# ranking and scaling data
def scale_data(data, max_poss):
    data_to_rank = data[data != 0]
    max_obs = np.max(data_to_rank)
    rank_data = rankdata(data_to_rank, method='average')
    rank_data = rank_data/(len(data_to_rank))
    rank_data = rank_data * max_obs/max_poss

    original_data_shape = data.shape
    data = data.astype(float).flatten()
    np.put(data, np.where(data != 0)[0], rank_data)
    data = data.reshape(original_data_shape)
    return data


def main(metric_path, gene_id_mapped_path, out_path, num_cpus):
    print('Loading data...', flush=True)
    # load metric data
    f = open(metric_path, 'r')
    metric_top1000 = []
    for line in f.readlines():
        line = line.strip().split('\t')
        metric_top1000.append(line)
    metric_top1000 = np.array(metric_top1000)
    f.close()

    # find the mapped gene ids that have metric associations listed
    unique_gene_codes_with_associations = set(metric_top1000[:, 0])

    f = open(gene_id_mapped_path, 'r')
    all_regions_gene_ids = []
    for line in f.readlines():
        region_gene_ids = []
        geneids = re.split('[;|\t]', line.strip())[3:]  # '\t|;|\|'
        for geneid in geneids:
            geneid = geneid.split('.')[0]  # removing everything after .
            if geneid in unique_gene_codes_with_associations:
                region_gene_ids.append(geneid)
        all_regions_gene_ids.append(region_gene_ids)
    f.close()

    print('Mapping genes to regions, genes, and scores...', flush=True)
    # making a dictionary to map genes to regions. Each gene will map to an
    # arrray of regions
    gene2region = map_genes_to_regions(all_regions_gene_ids)

    # making dictionaries to map a gene to similar genes and corresponding
    # scores
    gene2genes, gene2scores = map_genes_to_genes(metric_top1000)

    # split on cluster and find similar regions matrix based on STING metric
    # scores
    print('Finding similar regions...', flush=True)
    func = partial(create_sim_matrix, all_regions_gene_ids, metric_top1000,
                   gene2region, gene2genes, gene2scores)
    with closing(Pool(num_cpus)) as pool:
        top_1000_similar_inds = pool.map(func,
                                         np.arange(len(all_regions_gene_ids)))
    pool.join()
    top_1000_similar_inds = np.array(top_1000_similar_inds)
    # unsplit results from pool function
    top_1000_region_scores = top_1000_similar_inds[:, :1000].astype(int)
    top_1000_similar_regions = top_1000_similar_inds[:, 1000:2000].astype(int)
    top_1000_origin_gene_tracker = np.char.encode(
                                       top_1000_similar_inds[:, 2000:3000],
                                       'utf-8'
                                   )
    top_1000_destination_gene_tracker = np.char.encode(
                                            top_1000_similar_inds[:, 3000:],
                                            'utf-8'
                                        )
    del top_1000_similar_inds

    # Scaling similarity values and setting scores of score values=0 to 0
    # max_possible_val = 1000
    # min_nonzero_val = np.min(top_1000_region_scores[
    #                              top_1000_region_scores != 0
    #                          ])
    # top_1000_scaled_vals = ((top_1000_region_scores - min_nonzero_val) /
    #                         (max_possible_val - min_nonzero_val))
    # top_1000_scaled_vals[top_1000_scaled_vals<0] = 0
    max_poss = 1000
    top_1000_scaled_vals = scale_data(top_1000_region_scores, max_poss)

    # save top 100 of cosine similarities, similar regions, and percentiles
    #  for each region
    metric_f = h5py.File(out_path, 'w')
    metric_f.create_dataset('top_1000_scaled_vals', data=top_1000_scaled_vals)
    metric_f.create_dataset('top_1000_similar_regions',
                            data=top_1000_similar_regions)
    metric_f.create_dataset('top_1000_STRING_vals',
                            data=top_1000_region_scores)
    metric_f.create_dataset('top_1000_original_genes',
                            data=top_1000_origin_gene_tracker)
    metric_f.create_dataset('top_1000_destination_genes',
                            data=top_1000_destination_gene_tracker)
    metric_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
