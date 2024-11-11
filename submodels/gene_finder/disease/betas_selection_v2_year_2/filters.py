import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
plt.style.use("ggplot")
import math



def benjamini_hochberg(betas_matrix_average):
    betas_bin_average = np.mean(betas_matrix_average , axis = 0)
    betas_std = np.std(betas_matrix_average , axis = 0)
    z_score = betas_bin_average / betas_std #(np.std(betas_matrix_average , axis = 0) / np.sqrt(betas_matrix_average.shape[0]))
    p_value = stats.norm.sf(abs(z_score)) * 2
    p_value_sort = sorted(p_value)

    dict_beta_original_index = {}
    for i in range(len(p_value_sort)):
        for j in range(len(p_value)):
            if p_value_sort[i] == p_value[j]:
                dict_beta_original_index[i] = j
                break


    list_indexes = np.array(list(range(len(p_value_sort)))) + 1
    list_indexes = np.array(sorted(list_indexes))
    list_hochberg = 0.05 * (list_indexes / len(list_indexes))

    for i in range(len(list_indexes)):
        p_val = p_value_sort[i]
        hoch_val = list_hochberg[i]
        if p_val >= hoch_val:
            break

    list_significant_bins = list_indexes[:i]
    list_kept_index_betas = []
    for index in list_significant_bins:
        list_kept_index_betas.append(dict_beta_original_index[index - 1])
    #for i in list_significant_bins:
        #zeros_list = np.zeros(betas_matrix_average.shape[0])
        #betas_matrix_average[: , i - 1] = zeros_list
    #betas_matrix_average_new = np.delete(betas_matrix_average , np.s_[list_deleted_index_betas[-1] -1 : list_deleted_index_betas[0]] , axis = 1)

    return list_kept_index_betas

def bayesian_interval(betas_matrix_average):
    betas_matrix_new = np.empty(betas_matrix_average.shape)
    two_five_percent = math.ceil(betas_matrix_average.shape[0] * 0.025)
    ninety_five_percent = math.ceil(betas_matrix_average.shape[0] * 0.975)
    #get each column and order, then get percentiles and remove data
    for i in range(betas_matrix_average.shape[1]):
        bin_current = betas_matrix_average[:,i]
        bin_current = list(sorted(bin_current))
        betas_matrix_new[:,i] = bin_current

    list_delete_index = []
    betas_matrix_average_new = []
    list_keep_index = []
    for i in range(betas_matrix_new.shape[1]):
        if betas_matrix_new[two_five_percent , i] < 0 and betas_matrix_new[ninety_five_percent , i] > 0:
            list_delete_index.append(i)
        else:
            betas_matrix_average_new.append(betas_matrix_average[: , i])
            list_keep_index.append(i)


    #for index in list_delete_index:
     #   zeros_list = np.zeros(betas_matrix_average.shape[0])
     #   betas_matrix_average[:,index] = zeros_list



    return list_keep_index


def read_file(name_file):
    theta_vector_matrix = []
    for index , line in enumerate(open(name_file)):
        list_theta_vectors = []
        list_theta_strings = line.split("|")
        for theta in list_theta_strings:
            theta_vector = theta.strip().split(" ")
            if theta_vector != ['']:
                theta_vector = list(map(lambda x: float(x) , theta_vector))
                list_theta_vectors.append(theta_vector)
        theta_vector_matrix.append(list_theta_vectors)
    return np.array(theta_vector_matrix)


def average_theta_matrix(theta_matrix):
    (number_replica , number_theta_vectors_acum , number_theta_parameters) = theta_matrix.shape
    matrix_avg = np.zeros((number_theta_vectors_acum , number_theta_parameters))
    for i in range(number_replica):
        matrix_avg = matrix_avg + theta_matrix[i]
    matrix_avg = matrix_avg / number_replica
    return matrix_avg


def adjust_betas(betas_matrix , data_values , variance):
    #diag_var = np.diag(variance.reshape(len(variance)))
    var_inv = 1 / variance.reshape(len(variance))
    for i in range(len(var_inv)):
        if np.isinf(var_inv[i]):
            var_inv[i] = 0
    diag_var_inv = np.diag(var_inv)
    x_t_mul_diag_inv = np.matmul(data_values.T , diag_var_inv)
    a_matrix = np.matmul(x_t_mul_diag_inv , data_values)
    (eigval , eigvec) = np.linalg.eigh(a_matrix)
    rot = np.flip(eigvec , axis = 1)
    beta_star = np.matmul(rot.T , betas_matrix.T)
    for i in range(beta_star.shape[1]):
        for j in range(beta_star.shape[0]):
            if j > len(variance) - 1:
                beta_star[j , i] = 0

    beta_new = np.matmul(rot , beta_star)
    return(beta_new.T)


def read_data(file_matrix):
    list_matrix = []
    for index , line in enumerate(open(file_matrix)):
        data = line.strip().split(" ")
        data = list(map(lambda x:float(x) , data))
        list_matrix.append(data)
    return np.array(list_matrix)


def read_test_values(test_values):
    list_y_values = []
    for index , line in enumerate(open(test_values)):
        data = line.strip()
        list_y_values.append(float(data))
    return np.array(list_y_values).reshape((len(list_y_values) , 1))

