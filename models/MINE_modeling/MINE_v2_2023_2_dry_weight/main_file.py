import numpy as np
import monte_carlo as mc
import matrix_generation as mg
import utilities as util
import feature_selection as fs
import sys




#parameters of the simulation

dict_config_file = util.read_config_file(sys.argv[1])

name_simulation = dict_config_file["NAME_SIMULATION"]
file_matrix = dict_config_file["MATRIX_FILE"]
labels_file = dict_config_file["LABELS_FILE"]
number_replica = int(dict_config_file["NUMBER_REPLICA"])
number_equilibration_sweeps = int(dict_config_file["NUMBER_EQUILIBRATION_SWEEPS"])
number_decorrelation_sweeps = int(dict_config_file["NUMBER_DECORRELATION_SWEEPS"])
number_theta_vectors = int(dict_config_file["NUMBER_BETA_VECTORS"])
data_values = util.read_data(file_matrix)
(number_datapoints , number_betas) = data_values.shape
max_beta = 1.0
min_beta = 0.0
beta = [0.0] * number_betas
gamma = [0.0] * number_datapoints
step_width = np.ones(number_betas) * float(dict_config_file["STEP_WIDTH"])
delta_sw_adj_step = int(dict_config_file["STEP_WIDTH_ADJUSTMENT_INTERVAL"])
y_values = util.read_test_values(dict_config_file["PHENOTYPE_FILE"])
variance = util.read_test_values(dict_config_file["PHENOTYPE_VARIANCE_FILE"])
for i in range(len(variance)):
    if variance[i] == 0:
        variance[i] = 0.000000001
variance_rec = 1/variance
initial_beta = np.random.uniform(0 , 1 , (number_betas , 1)) #* (max_beta - min_beta) + min_beta
initial_gamma = np.random.uniform(0 , 1 , (number_datapoints , 1))
number_bins = int(dict_config_file["NUMBER_BINS"])
binwidth = (2 * max_beta) / number_bins
s_cut_evals = float(dict_config_file["THRESHOLD_EVALS"])
number_new_experiments = int(dict_config_file["NUMBER_NEW_EXPERIMENTS"])
number_max_elements_per_partition = int(dict_config_file["SUB_OPTIMAL_SEARCH_TUPLE"])
number_parallel_process_mine = int(dict_config_file["NUMBER_PARALLEL_PROCESSES_MINE"])
option_checkpoint = int(dict_config_file["OPTION_CHECKPOINT_START"])
time_checkpoint_suboptimal = int(dict_config_file["TIME_CHECKPOINT_SUBOPTIMAL"])
number_steps_monte_carlo_mine = int(dict_config_file["NUMBER_STEPS_MONTE_CARLO_MINE"])
criteria_option = int(dict_config_file["MINE_CRITERIA_OPTION"]) # 1 for correlation criteria, 2 for covariance criteria
option_wf = int(dict_config_file["OPTION_WORKFLOW"]) # 1 for simulation, 2 for MINE base, 3 for MINE n choose 3, 4 for MINE greedy, 5 for MINE nc3_2, 6 for MINE monte carlo



if option_wf == 1:
    (theta_vector_matrix , gamma_vector_matrix) = mc.mc_simulation(number_replica , initial_beta , initial_gamma , data_values , y_values , variance_rec , number_equilibration_sweeps , number_betas , step_width , max_beta , min_beta , number_theta_vectors , number_decorrelation_sweeps , delta_sw_adj_step , name_simulation)
else:
    beta_vector_tensor = util.read_file(name_simulation + "/beta")
    gamma_vector_tensor = util.read_file(name_simulation + "/gamma")
    beta_vector_matrix_avg = util.average_theta_matrix(beta_vector_tensor)
    gamma_vector_matrix_avg = util.average_theta_matrix(gamma_vector_tensor)
    beta_matrix_star = fs.adjust_betas(beta_vector_matrix_avg, data_values, variance)
    labels = mg.get_labels(labels_file)
    data_values_range = mg.create_grid(data_values, number_max_elements_per_partition)
    if option_wf == 2:
        completion = mg.mine_criteria_suboptimal(beta_matrix_star , gamma_vector_matrix_avg , data_values_range , labels , data_values , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation , number_parallel_process_mine , option_checkpoint , time_checkpoint_suboptimal)
    elif option_wf == 3:
        file_accessions = mg.get_experiment_accessions_suboptimal(name_simulation, number_new_experiments , criteria_option)
        mg.get_mine_score_final_accessions(file_accessions, labels, data_values, beta_matrix_star, gamma_vector_matrix_avg, number_datapoints, s_cut_evals, number_new_experiments, name_simulation, "labels_score_Nc3")
    elif option_wf == 4:
        file_accessions = mg.mine_criteria_greedy(beta_matrix_star , gamma_vector_matrix_avg , labels , data_values , number_datapoints , s_cut_evals , number_new_experiments , name_simulation , criteria_option)
        mg.get_mine_score_final_accessions(file_accessions, labels, data_values, beta_matrix_star, gamma_vector_matrix_avg, number_datapoints, s_cut_evals, number_new_experiments, name_simulation, "labels_score_greedy")
    elif option_wf == 5:
        (labels_nc3 , indexes_nc3) = mg.get_labels_suboptimal(name_simulation , criteria_option)
        file_accessions = mg.mine_criteria_nc3_2(beta_vector_matrix_avg , gamma_vector_matrix_avg , labels , labels_nc3 , indexes_nc3 , data_values , number_datapoints , s_cut_evals , number_max_elements_per_partition , name_simulation , number_new_experiments , criteria_option)
        mg.get_mine_score_final_accessions(file_accessions, labels, data_values, beta_matrix_star, gamma_vector_matrix_avg, number_datapoints, s_cut_evals, number_new_experiments, name_simulation, "labels_score_Nc3_Nc2")
    elif option_wf == 6:
        mg.mine_criteria_monte_carlo(beta_vector_matrix_avg, gamma_vector_matrix_avg, labels, data_values, number_datapoints, s_cut_evals, number_new_experiments, name_simulation, criteria_option, number_steps_monte_carlo_mine , "labels_score_monte_carlo")
