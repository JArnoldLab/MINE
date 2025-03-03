import filters
import numpy as np
import genes_retrieval as gt


data_values = filters.read_data("matrix_dry_weight_log20_v3")
variance = filters.get_initial_sigma_square("variance")
betas_matrix_mc = filters.read_file("beta")
betas_matrix_mc_average = filters.average_theta_matrix(betas_matrix_mc)
betas_matrix_star = betas_matrix_mc_average#filters.adjust_betas(betas_matrix_mc_average , data_values , variance)
phenotype = "dry_weight"
name_simulation = "mixed"

index_kept_bayesian_interval = filters.bayesian_interval(betas_matrix_star)

index_kept_benjamini_hochberg = filters.benjamini_hochberg(betas_matrix_star)

set_kept_bayesian_interval = set(index_kept_bayesian_interval)
set_kept_benjamini_hochberg = set(index_kept_benjamini_hochberg)
index_intersection = set_kept_bayesian_interval.intersection(set_kept_benjamini_hochberg)
index_union_unique = set_kept_bayesian_interval.union(set_kept_benjamini_hochberg)
bins_remaining_index = list(sorted(index_intersection))
gt.find_genes(bins_remaining_index , phenotype , name_simulation)
