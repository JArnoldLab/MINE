import numpy as np
import operator
import itertools
from scipy import linalg
import time
import random


def calculate_expected(g_function):
    g_function_sum = np.sum(g_function)
    return g_function_sum / len(g_function)


def calculate_expected_two_values(g_function_1 , g_function_2):
    g_function_mult = g_function_1 * g_function_2
    g_function_sum = np.sum(g_function_mult)
    return g_function_sum / len(g_function_mult)


def d_matrix_creation(beta_matrix , gamma_matrix , data_values , data_vectors , number_new_experiments):
    d_matrix = np.empty((number_new_experiments , number_new_experiments))
    counter_u = 0
    g_function_list = []
    for i in range(number_new_experiments):
        g_function = np.dot(beta_matrix , data_vectors[data_values[i]])
        g_function_list.append(g_function)

    for i in range(number_new_experiments):
        for j in range(number_new_experiments):
            experiment = i
            experiment_prime = j
            g_function_i = g_function_list[experiment]
            g_function_j = g_function_list[experiment_prime]

            expected_two_values = calculate_expected_two_values(g_function_i , g_function_j)
            expected_one = calculate_expected(g_function_i)
            expected_two = calculate_expected(g_function_j)
            value = expected_two_values - (expected_one * expected_two)
            d_matrix[i][j] = value
    return d_matrix


def e_matrix_creation(d_matrix):
    (rows_number , cols_number) = d_matrix.shape
    e_matrix = np.empty((rows_number , cols_number))
    for i in range(rows_number):
        for j in range(cols_number):
            if d_matrix[i][i] != 0 and d_matrix[j][j] != 0:
                e_matrix[i][j] = d_matrix[i][j] / (np.sqrt(d_matrix[i][i]) * np.sqrt(d_matrix[j][j]))
            else:
                e_matrix[i][j] = 0

    return e_matrix


def regularize_evals(evals , s_cut):
    min_value = evals[0] / s_cut
    for i in range(len(evals)):
        evals[i] = max(evals[i] , min_value)
    return evals


def determinant_squared_diag(eivals):
    log_eivals = np.log10(eivals)
    determinant = np.sum(log_eivals)
    return determinant

def print_d_matrix_e_format(d_matrix):
    file_output = open("d_matrix.csv", "w")
    file_output.write("index_row,index_col,element_matrix\n")
    for j in range(len(d_matrix[0])):
        for i in range(len(d_matrix)):
            file_output.write(str(i) + "," + str(j) + "," + str(d_matrix[i][j]) + "\n")
    file_output.close()

def print_evals_e_format(evals):
    evals = sorted(evals , reverse = True)
    file_output = open("evals.csv" , "w")
    file_output.write("index,elemennt_array\n")
    for i in range(len(evals)):
        file_output.write(str(i) + "," + str(evals[i]) + "\n")
    file_output.close()

def matrix_verification(d_matrix , evals , evecs):
    list_difference = []
    max_eigval = max(abs(evals))
    for i in range(len(evecs[0])):
        matrix_part = np.matmul(d_matrix , evecs[:,i])
        evecs_evals_part = evals[i] * evecs[:,i]
        difference = matrix_part - evecs_evals_part
        rel_err = np.linalg.norm(difference) / (max_eigval * np.linalg.norm(evecs[:,i]))
        list_difference.append(rel_err)
    return list_difference

def read_bernd_file_monte_carlo(bernd_file_montecarlo):
    list_theta = []
    list_theta_vector = []
    for index , line in enumerate(open(bernd_file_montecarlo)):
        list_theta = []
        if index != 0:
            data = line.split(",")
            list_theta.append(float(data[3].strip()))
            list_theta.append(float(data[4].strip()))
            list_theta.append(float(data[5].strip()))
            list_theta_vector.append(list_theta)
    return np.array(list_theta_vector)

def read_bernd_file_u_vectors(bernd_file_u_vectors):
    list_experiments = []
    list_experiment_1 = []
    list_experiment_2 = []
    bernd_results = {}
    for index , line in enumerate(open(bernd_file_u_vectors)):
        list_experiment_1 = []
        list_experiment_2 = []
        if index != 0:
            data = line.split(",")
            list_experiment_1.append(float(data[3].strip()))
            list_experiment_1.append(float(data[4].strip()))
            list_experiment_1.append(float(data[5].strip()))
            list_experiment_2.append(float(data[6].strip()))
            list_experiment_2.append(float(data[7].strip()))
            list_experiment_2.append(float(data[8].strip()))
            list_experiments.append((list_experiment_1 , list_experiment_2))
            bernd_results[str(list_experiment_1[0]) + "," + str(list_experiment_1[1]) + "," + str(list_experiment_1[2]) + "," + str(list_experiment_2[0]) + "," + str(list_experiment_2[1]) + "," + str(list_experiment_2[2])] = float(data[1].strip())

    return (list_experiments , bernd_results)



def mine_criteria_core_suboptimal(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation):


    d_matrix_u_set = d_matrix_creation(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , data_vectors , number_new_experiments)
    e_matrix_u_set = e_matrix_creation(d_matrix_u_set)
    evals_d_matrix = sorted(linalg.eigh(d_matrix_u_set , eigvals_only = True) , reverse = True)
    evals_e_matrix = sorted(linalg.eigh(e_matrix_u_set , eigvals_only = True) , reverse = True)

    evals_d_matrix = regularize_evals(evals_d_matrix, s_cut_evals)
    evals_e_matrix = regularize_evals(evals_e_matrix, s_cut_evals)
    result_d_matrix = determinant_squared_diag(evals_d_matrix)
    result_e_matrix = determinant_squared_diag(evals_e_matrix)


    file_output = open(name_simulation + "/labels_results" , "a")

    for i in range(len(data_values_range)):
        file_output.write(labels[data_values_range[i]] + " ")
    file_output.write("|" + str(result_d_matrix) + "|" + str(result_e_matrix) + "\n")
    file_output.close()




def get_labels(labels_file):
    list_labels = []
    for index , line in enumerate(open(labels_file)):
        list_labels.append(line.strip())

    return (list_labels)


def create_grid(data_vectors , number_max_partitions):
    list_indexes = list(range(len(data_vectors)))
    list_combinations = itertools.combinations(list_indexes , number_max_partitions)

    return list_combinations

def read_grid(file_grid):
    list_data_values = []
    for index , line in enumerate(open(file_grid)):
        data = line.strip().split("|")
        list_data_values.append(data)
    return list_data_values


def get_check_point(data_value_tuple , name_simulation):
    file_checkpoint = open(name_simulation + "/checkpoint" , "w")
    for item in data_value_tuple:
        file_checkpoint.write(str(item) + ",")
    file_checkpoint.close()

def read_checkpoint(name_simulation):
    list_elements = []
    for index , line in enumerate(open(name_simulation + "/checkpoint")):
        data = line.split(",")
        list_elements = data[:len(data) - 1]
    list_elements = list(map(lambda x: int(x) , list_elements))
    return(tuple(list_elements))



def mine_criteria_suboptimal(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , labels , data_vectors , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation , number_parallel_process_mine , option_checkpoint , time_checkpoint_suboptimal):

    if option_checkpoint == 1:
        data_value_tuple_check = read_checkpoint(name_simulation)
        time_ini = time.perf_counter()
        flag_checkpoint = 0
        for data_value_tuple in data_values_range:
            if data_value_tuple == data_value_tuple_check:
                flag_checkpoint = 1
            if flag_checkpoint == 1:
                mine_criteria_core_suboptimal(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_value_tuple , labels , data_vectors , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation)
                time_fin = time.perf_counter()
                if time_fin - time_ini >= time_checkpoint_suboptimal:
                    get_check_point(data_value_tuple, name_simulation)
                    return -1
    else:
        file_output_labels = open(name_simulation + "/labels_results", "w")
        file_output_labels.write("Accessions|log(det) D-matrix|log(det) E-matrix\n")
        file_output_labels.close()
        time_ini = time.perf_counter()
        for data_value_tuple in data_values_range:
            mine_criteria_core_suboptimal(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_value_tuple , labels , data_vectors , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation)
            time_fin = time.perf_counter()
            if time_fin - time_ini >= time_checkpoint_suboptimal:
                get_check_point(data_value_tuple , name_simulation)
                return -1


    return 1



def get_experiment_accessions_suboptimal(name_simulation , number_new_experiments , criteria_option):
    list_data_d_matrix = []
    list_data_e_matrix = []
    for index , line in enumerate(open(name_simulation + "/labels_results")):
        if index != 0:
            data = line.split("|")
            list_data_d_matrix.append((data[0] , float(data[1])))
            list_data_e_matrix.append((data[0] , float(data[2].strip())))
    list_data_d_matrix = sorted(list_data_d_matrix , key = operator.itemgetter(1) , reverse = True)
    list_data_e_matrix = sorted(list_data_e_matrix , key = operator.itemgetter(1) , reverse = True)

    cont_elements = 0
    list_repeated = []
    file_output = open(name_simulation + "/final_accessions_d_matrix" , "w")
    file_output_top = open(name_simulation + "/top_results_d" , "w")
    for i in range(len(list_data_d_matrix)):
        data = list_data_d_matrix[i][0].strip().split(" ")
        file_output_top.write(str(list_data_d_matrix[i][0]) + "|" + str(list_data_d_matrix[i][1]) + "\n")
        for j in range(len(data)):
            if data[j] not in list_repeated:
                file_output.write(data[j] + "\n")
                cont_elements = cont_elements + 1
                list_repeated.append(data[j])
                if cont_elements >= number_new_experiments:
                    break
                list_repeated.append(data[j])
        if cont_elements >= number_new_experiments:
            break

    file_output.close()
    file_output_top.close()

    cont_elements = 0
    list_repeated = []
    file_output = open(name_simulation + "/final_accessions_e_matrix" , "w")
    file_output_top = open(name_simulation + "/top_results_e" , "w")
    for i in range(len(list_data_e_matrix)):
        data = list_data_e_matrix[i][0].strip().split(" ")
        file_output_top.write(str(list_data_e_matrix[i][0]) + "|" + str(list_data_e_matrix[i][1]) + "\n")
        for j in range(len(data)):
            if data[j] not in list_repeated:
                file_output.write(data[j] + "\n")
                cont_elements = cont_elements + 1
                list_repeated.append(data[j])
                if cont_elements >= number_new_experiments:
                    break
        if cont_elements >= number_new_experiments:
            break

    file_output.close()
    file_output_top.close()

    if criteria_option == 1:
        file_name = name_simulation + "/final_accessions_e_matrix"
    else:
        file_name = name_simulation + "/final_accessions_d_matrix"

    return file_name



def mine_criteria_core_greedy(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation , file_prefix):
    list_values_det = []
    list_values_det_e = []
    #for i in range(len(data_values_range)):
    d_matrix_u_set = d_matrix_creation(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , data_vectors , number_new_experiments)
    e_matrix_u_set = e_matrix_creation(d_matrix_u_set)
    evals_d_matrix = sorted(linalg.eigh(d_matrix_u_set , eigvals_only = True) , reverse = True)
    evals_e_matrix = sorted(linalg.eigh(e_matrix_u_set , eigvals_only = True) , reverse = True)

    evals_d_matrix = regularize_evals(evals_d_matrix, s_cut_evals)
    evals_e_matrix = regularize_evals(evals_e_matrix, s_cut_evals)
    result_d_matrix = determinant_squared_diag(evals_d_matrix)
    result_e_matrix = determinant_squared_diag(evals_e_matrix)
        #list_values_det.append((result_d_matrix , i))
        #list_values_det_e.append((result_e_matrix , i))

    #list_values_det = sorted(list_values_det , key = operator.itemgetter(0) , reverse = True)
    #list_values_det_e = sorted(list_values_det_e , key = operator.itemgetter(0) , reverse = True)
    #print(data_values_range)
    #print("\n")
    #print("New subset")
    file_output = open(name_simulation + "/labels_results_" + file_prefix , "a")

    for i in range(len(data_values_range)):
        file_output.write(labels[data_values_range[i]] + " ")
    file_output.write("|" + str(result_d_matrix) + "|" + str(result_e_matrix) + "\n")

    file_output.close()


def mine_criteria_greedy(beta_vector_matrix_avg , gamma_vector_matrix_avg , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation , criteria_option):

    top_tuple = get_experiment_accessions_greedy(name_simulation , labels , "main" , criteria_option)
    counter_items = len(top_tuple)
    while counter_items < number_new_experiments:
        file_output_labels_greedy = open(name_simulation + "/labels_results_sub", "w")
        file_output_labels_greedy.write("Accessions|log(det) D-matrix|log(det) E-matrix\n")
        file_output_labels_greedy.close()
        for i in range(len(data_vectors)):
            if i not in top_tuple:
                top_tuple.extend(i)
                mine_criteria_core_greedy(beta_vector_matrix_avg , gamma_vector_matrix_avg , top_tuple , labels , data_vectors , number_datapoints , s_cut_evals, len(top_tuple) , name_simulation ,"sub")
                del top_tuple[-1]


        top_tuple = get_experiment_accessions_greedy(name_simulation , "sub" , labels)
        counter_items = counter_items + 1

    if criteria_option == 1:
        file_name = name_simulation + "/final_accessions_e_matrix_greedy"
        file_output = open(file_name , "w")
    else:
        file_name = name_simulation + "/final_accessions_d_matrix_greedy"
        file_output = open(file_name , "w")
    for i in range(len(top_tuple)):
        file_output.write(labels[i] + "\n")
    file_output.close()

    return file_name




def get_experiment_accessions_greedy(name_simulation , labels , prefix , criteria_option):
    list_data_d_matrix = []
    list_data_e_matrix = []
    list_data_indexes_e = []
    list_data_indexes_d = []

    if prefix == "main":
        prefix_total = "/labels_results"
    else:
        prefix_total = "/labels_results_" + prefix

    for index , line in enumerate(open(name_simulation + prefix_total)):
        if index != 0:
            data = line.split("|")
            list_data_d_matrix.append((data[0] , float(data[1])))
            list_data_e_matrix.append((data[0] , float(data[2].strip())))


    list_data_d_matrix = sorted(list_data_d_matrix , key = operator.itemgetter(1) , reverse = True)
    list_data_e_matrix = sorted(list_data_e_matrix , key = operator.itemgetter(1) , reverse = True)

    list_labels_e = list_data_e_matrix[0].split(",")
    list_labels_d = list_data_d_matrix[0].split(",")

    for item in list_labels_e:
        index = labels.index(item)
        list_data_indexes_e.append(index)

    for item in list_labels_d:
        index = labels.index(item)
        list_data_indexes_d.append(index)

    if criteria_option == 1:
        return list_data_indexes_e
    else:
        return list_data_indexes_d



def mine_criteria_core_nc3_2(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation):
    list_values_det = []
    list_values_det_e = []
    #for i in range(len(data_values_range)):
    d_matrix_u_set = d_matrix_creation(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , data_vectors , number_new_experiments)
    e_matrix_u_set = e_matrix_creation(d_matrix_u_set)
    evals_d_matrix = sorted(linalg.eigh(d_matrix_u_set , eigvals_only = True) , reverse = True)
    evals_e_matrix = sorted(linalg.eigh(e_matrix_u_set , eigvals_only = True) , reverse = True)

    evals_d_matrix = regularize_evals(evals_d_matrix, s_cut_evals)
    evals_e_matrix = regularize_evals(evals_e_matrix, s_cut_evals)
    result_d_matrix = determinant_squared_diag(evals_d_matrix)
    result_e_matrix = determinant_squared_diag(evals_e_matrix)
        #list_values_det.append((result_d_matrix , i))
        #list_values_det_e.append((result_e_matrix , i))

    #list_values_det = sorted(list_values_det , key = operator.itemgetter(0) , reverse = True)
    #list_values_det_e = sorted(list_values_det_e , key = operator.itemgetter(0) , reverse = True)
    #print(data_values_range)
    #print("\n")
    #print("New subset")
    file_output = open(name_simulation + "/labels_results" , "a")

    for i in range(len(data_values_range)):
        file_output.write(labels[data_values_range[i]] + " ")
    file_output.write("|" + str(result_d_matrix) + "|" + str(result_e_matrix) + "\n")
    file_output.close()


def mine_criteria_nc3_2(beta_vector_matrix_avg , gamma_vector_matrix_avg , labels , labels_nc3 , data_vectors , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation , number_new_experiments , criteria_option):
    list_accumulate = []
    list_accumulate_index = []
    labels_current = labels.copy()
    data_vectors_current = data_vectors.copy()
    list_accumulate.extend(labels_nc3)
    (data_vectors_current , labels_current) = remove_labels_data(data_vectors_current , labels_current , list_accumulate)
    i = 3
    while i < number_new_experiments:
        data_values_range = create_grid(data_vectors_current , number_max_elements_per_partition)
        number_elements = number_max_elements_per_partition

        file_output_labels = open(name_simulation + "/labels_results", "w")
        file_output_labels.write("Accessions|log(det) D-matrix|log(det) E-matrix\n")
        file_output_labels.close()

        for data_values in data_values_range:
            mine_criteria_core_nc3_2(beta_vector_matrix_avg, gamma_vector_matrix_avg, data_values, labels, data_vectors, number_datapoints, s_cut_evals, number_elements, name_simulation)

        (labels_d , labels_e , log_det_d , log_det_e) = get_experiment_accessions_nc3_2(name_simulation)
        if criteria_option == 1:
            list_accumulate.extend(labels_e)
            (data_vectors_current , labels_current) = remove_labels_data(data_vectors_current , labels_current , labels_e)
        else:
            list_accumulate.extend(labels_d)
            (data_vectors_current , labels_current) = remove_labels_data(data_vectors_current , labels_current , labels_d)
        i = i + number_elements


    if criteria_option == 1:
        file_name = name_simulation + "/final_accessions_e_matrix"
        file_output_final_labels = open(file_name , "w")
    elif criteria_option == 2:
        file_name = name_simulation + "/final_accessions_d_matrix"
        file_output_final_labels = open(file_name , "w")

    for i in range(len(list_accumulate)):
        file_output_final_labels.write(list_accumulate[i] + "\n")
    file_output_final_labels.close()

    return file_name




def remove_labels_data(data_vectors , labels , labels_remove):
    list_indexes = []
    for element in labels_remove:
        index_element = labels.index(element)
        del labels[index_element]
        data_vectors = np.delete(data_vectors , index_element , 0)

    return (data_vectors , labels)

def get_experiment_accessions_nc3_2(name_simulation):
    list_data_d_matrix = []
    list_data_e_matrix = []
    for index , line in enumerate(open(name_simulation + "/labels_results")):
        if index != 0:
            data = line.split("|")
            list_data_d_matrix.append((data[0] , float(data[1])))
            list_data_e_matrix.append((data[0] , float(data[2].strip())))
    list_data_d_matrix = sorted(list_data_d_matrix , key = operator.itemgetter(1) , reverse = True)
    list_data_e_matrix = sorted(list_data_e_matrix , key = operator.itemgetter(1) , reverse = True)

    return(list_data_d_matrix[0][0].strip().split(" ") , list_data_e_matrix[0][0].strip().split(" ") , list_data_d_matrix[0][1] , list_data_e_matrix[0][1])



def mine_criteria_core_monte_carlo(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation):
    list_values_det = []
    list_values_det_e = []
    #for i in range(len(data_values_range)):
    d_matrix_u_set = d_matrix_creation(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_values_range , data_vectors , number_new_experiments)
    e_matrix_u_set = e_matrix_creation(d_matrix_u_set)
    evals_d_matrix = sorted(linalg.eigh(d_matrix_u_set , eigvals_only = True) , reverse = True)
    evals_e_matrix = sorted(linalg.eigh(e_matrix_u_set , eigvals_only = True) , reverse = True)

    evals_d_matrix = regularize_evals(evals_d_matrix, s_cut_evals)
    evals_e_matrix = regularize_evals(evals_e_matrix, s_cut_evals)
    result_d_matrix = determinant_squared_diag(evals_d_matrix)
    result_e_matrix = determinant_squared_diag(evals_e_matrix)
        #list_values_det.append((result_d_matrix , i))
        #list_values_det_e.append((result_e_matrix , i))

    #list_values_det = sorted(list_values_det , key = operator.itemgetter(0) , reverse = True)
    #list_values_det_e = sorted(list_values_det_e , key = operator.itemgetter(0) , reverse = True)
    return(result_d_matrix , result_e_matrix)



def mine_criteria_monte_carlo(beta_vector_matrix_avg , gamma_vector_matrix_avg , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation , criteria_option , number_steps_annealing , prefix):
    file_output_labels = open(name_simulation + "/" + prefix + "_labels_results", "w")
    file_output_chi = open(name_simulation + "/chi_sq" , "w")
    file_output_labels.write("Accessions|log(det) D-matrix|log(det) E-matrix\n")
    file_output_chi.write("log(det) D-matrix|log(det) E-matrix|step\n")


    data = select_random_accessions(data_vectors , number_new_experiments)
    (log_det_d , log_det_e) = mine_criteria_core_monte_carlo(beta_vector_matrix_avg , gamma_vector_matrix_avg , data , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation)
    for steps in range(number_steps_annealing):
        data_new = data.copy()
        index_replacement = random.randint(0 , len(data) - 1)
        index_new = random.randint(0 , len(data_vectors) - 1)
        data_new[index_replacement] = index_new
        (log_det_d_new , log_det_e_new) = mine_criteria_core_monte_carlo(beta_vector_matrix_avg , gamma_vector_matrix_avg , data_new , labels , data_vectors , number_datapoints , s_cut_evals , number_new_experiments , name_simulation)
        if criteria_option == 1:
            delta_log_det = -(log_det_e_new) - -(log_det_e)
        else:
            delta_log_det = -(log_det_d_new) - -(log_det_d)

        if delta_log_det < 0:
            ratio = 2.0
        else:
            ratio = np.exp(-delta_log_det)

        if ratio >= 1:
            data = data_new
            log_det_e = log_det_e_new
            log_det_d = log_det_d_new
        else:
            v = np.random.uniform(0 , 1)
            if v < ratio:
                data = data_new
                log_det_e = log_det_e_new
                log_det_d = log_det_d_new

        file_output_chi.write(str(log_det_d) + "|" + str(log_det_e) + "|" + str(steps) + "\n")


    for i in range(len(data)):
        if i != len(data) - 1:
            file_output_labels.write(labels[data[i]] + " ")
        else:
            file_output_labels.write(labels[data[i]])
    file_output_labels.write("|" + str(log_det_d) + "|" + str(log_det_e) + "\n")
    file_output_labels.close()
    file_output_chi.close()


def select_random_accessions(data_vectors , number_new_experiments):
    list_elements = list(range(len(data_vectors)))
    list_selection = random.sample(list_elements , number_new_experiments)
    return list_selection



def get_labels_suboptimal(name_simulation , criteria_option):
    list_data_d_matrix = []
    list_data_e_matrix = []
    for index, line in enumerate(open(name_simulation + "/labels_results")):
        if index != 0:
            data = line.split("|")
            list_data_d_matrix.append((data[0], float(data[1])))
            list_data_e_matrix.append((data[0], float(data[2].strip())))
    list_data_d_matrix = sorted(list_data_d_matrix, key=operator.itemgetter(1), reverse=True)
    list_data_e_matrix = sorted(list_data_e_matrix, key=operator.itemgetter(1), reverse=True)

    if criteria_option == 1:
        labels_return = list_data_e_matrix[0][0].split(" ")
    else:
        labels_return = list_data_d_matrix[0][0].split(" ")

    return labels_return


def get_mine_score_final_accessions(file_accessions , labels , data_vectors , beta_matrix , gamma_matrix , number_datapoints , s_cut_evals , number_new_experiments , name_simulation , prefix):
    list_indexes_labels = []
    for index , line in enumerate(open(file_accessions)):
        index = labels.index(line.split())
        list_indexes_labels.append(index)

    mine_criteria_core_suboptimal(beta_matrix, gamma_matrix, list_indexes_labels, labels, data_vectors, number_datapoints, s_cut_evals, number_new_experiments, name_simulation + "/" + prefix)



