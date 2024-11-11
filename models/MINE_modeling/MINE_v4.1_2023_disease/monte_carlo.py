import numpy as np
import matplotlib.pyplot as plt
import random
import time
import utilities as util
from scipy import linalg

plt.style.use('ggplot')


def hamiltonian_func(y_values , data_values , beta , V_matrix_inverse , log_determinant_V):
    y_hat = np.dot(data_values , beta)
    y_diff = y_values - y_hat
    first_term = 0.5 * np.sum(np.square(y_diff.reshape(len(y_diff))) * V_matrix_inverse)
    second_term = ((len(y_values) * 0.5) * np.log(2 * np.pi)) + (0.5 * log_determinant_V)
    hamiltonian = first_term + second_term
    return (hamiltonian , y_hat , y_diff)


def hamiltonian_func_fast(y_diff , V_matrix_inverse , log_determinant_V):
    first_term = 0.5 * np.sum(np.square(y_diff.reshape(len(y_diff))) * V_matrix_inverse)
    second_term = ((len(y_diff) * 0.5) * np.log(2 * np.pi)) + (0.5 * log_determinant_V)
    hamiltonian = first_term + second_term
    return (hamiltonian , first_term , second_term)





def ensemble(delta_hamiltonian):
    return np.exp(-delta_hamiltonian)


def mc_update(beta , gamma , sigma_square , rho , hamiltonian , y_values , y_hat , y_diff , variance_prediction , delta_variance_prediction , number_parameters , delta_beta , delta_rho , data_values , V_matrix_inverse , log_determinant_V , variance_info , max_beta , min_beta , theta_selected , dict_variance_order , dict_variance_block , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff):
    acceptance = 0
    if theta_selected < len(beta):
        beta_new = beta.copy()
        gamma_new = gamma.copy()
        beta_new[theta_selected] = beta_new[theta_selected] + delta_beta  # in the general case I would alter on the chosen variable
        acceptance_acum_vector_fixed_eff[theta_selected] = acceptance_acum_vector_fixed_eff[theta_selected] + 1
        delta_y_hat = (data_values[:,theta_selected] * delta_beta).reshape(y_values.shape)
        y_hat_new = y_hat + delta_y_hat
        y_diff_new = y_values - y_hat_new
        variance_prediction_new = np.sum(np.square(y_diff_new)) / (len(y_diff_new) - 1)
        delta_variance_prediction_new = variance_prediction_new - variance_prediction
        (hamiltonian_new , first_term_new , second_term_new) = hamiltonian_func_fast(y_diff_new , V_matrix_inverse , log_determinant_V)
    else:
        sigma_square_new = sigma_square.copy()
        rho_new = rho.copy()
        rho_selected = theta_selected - len(beta)
        block_selected = dict_variance_block[rho_selected]
        rho_sigma_change = dict_variance_order[block_selected]
        acceptance_acum_vector_random_eff[rho_selected] = acceptance_acum_vector_random_eff[rho_selected] + 1
        for rho_sigma_change_selected in rho_sigma_change:
            rho_new[rho_sigma_change_selected] = rho_new[rho_sigma_change_selected] + delta_rho
            sigma_square_new[rho_sigma_change_selected] = np.exp(rho_new[rho_sigma_change_selected])
        (V_matrix_inverse_new , log_determinant_V_new) = update_V_matrix(data_values , sigma_square_new , rho_sigma_change , V_matrix_inverse , variance_prediction)
        (hamiltonian_new , first_term_new , second_term_new) = hamiltonian_func_fast(y_diff , V_matrix_inverse_new , log_determinant_V_new)



    delta_hamiltonian = hamiltonian_new - hamiltonian

    if delta_hamiltonian < 0:
        ratio = 2.0
    else:
        ratio = ensemble(delta_hamiltonian)

    if ratio >= 1:
        acceptance = 1
        if theta_selected < len(beta):
            acceptance_rate_vector_fixed_eff[theta_selected] = acceptance_rate_vector_fixed_eff[theta_selected] + 1
        else:
            acceptance_rate_vector_random_eff[rho_selected] = acceptance_rate_vector_random_eff[rho_selected] + 1
    else:
        v = np.random.uniform(0, 1)
        if v < ratio:
            acceptance = 1
            if theta_selected < len(beta):
                acceptance_rate_vector_fixed_eff[theta_selected] = acceptance_rate_vector_fixed_eff[theta_selected] + 1
            else:
                acceptance_rate_vector_random_eff[rho_selected] = acceptance_rate_vector_random_eff[rho_selected] + 1




    if acceptance == 1:
        if theta_selected < len(beta):
            return (beta_new , gamma_new , sigma_square , rho , hamiltonian_new , y_hat_new , V_matrix_inverse , log_determinant_V , y_diff_new , variance_prediction_new , delta_variance_prediction_new , first_term_new , second_term_new , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff)
        else:
            return (beta , gamma , sigma_square_new , rho_new , hamiltonian_new , y_hat , V_matrix_inverse_new , log_determinant_V_new , y_diff , variance_prediction , delta_variance_prediction , first_term_new , second_term_new , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff)
    else:
        return (beta , gamma , sigma_square , rho , hamiltonian , y_hat , V_matrix_inverse , log_determinant_V , y_diff , variance_prediction , delta_variance_prediction , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff)



def step_width_adjustment(acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_rate_vector_fixed_eff_old , acceptance_rate_vector_random_eff_old , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff , sw_rate_fixed_eff , sw_rate_random_eff , step_width_fixed_eff , step_width_random_eff , sw_adj_factor_fixed_eff , sw_adj_factor_random_eff , sw_low_accep_rate , sw_high_accep_rate , delta_sw_adj_step , equil_sweeps , step_width_fixed_eff_list , step_width_random_eff_list):
    acceptance_rate_vector_fixed_eff = acceptance_rate_vector_fixed_eff / acceptance_acum_vector_fixed_eff
    acceptance_rate_vector_random_eff = acceptance_rate_vector_random_eff / acceptance_acum_vector_random_eff
    sw_min_fixed_eff = sw_rate_fixed_eff * np.mean(np.array(step_width_fixed_eff_list) , axis = 0)
    sw_max_fixed_eff = (1 / sw_rate_fixed_eff) * np.mean(np.array(step_width_fixed_eff_list) , axis = 0)
    sw_min_random_eff = sw_rate_random_eff * np.mean(np.array(step_width_random_eff_list) , axis = 0)
    sw_max_random_eff = (1 / sw_rate_random_eff) * np.mean(np.array(step_width_random_eff_list) , axis = 0)
    for i in range(len(sw_adj_factor_fixed_eff)):
        if equil_sweeps == delta_sw_adj_step:
            if (acceptance_rate_vector_fixed_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_fixed_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_fixed_eff[i] < sw_low_accep_rate) or (acceptance_rate_vector_fixed_eff[i] < sw_low_accep_rate):
                sw_adj_factor_fixed_eff[i] = 2 / 3
            elif (acceptance_rate_vector_fixed_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_fixed_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_fixed_eff[i] > sw_high_accep_rate) or (acceptance_rate_vector_fixed_eff[i] > sw_high_accep_rate):
                sw_adj_factor_fixed_eff[i] = 3 / 2
        else:
            if (acceptance_rate_vector_fixed_eff[i] <= sw_high_accep_rate and acceptance_rate_vector_fixed_eff[i] >= sw_low_accep_rate) or (step_width_fixed_eff[i] < sw_min_fixed_eff[i]) or (step_width_fixed_eff[i] > sw_max_fixed_eff[i]):
                sw_adj_factor_fixed_eff[i] = 1
            elif (acceptance_rate_vector_fixed_eff_old[i] < sw_low_accep_rate and acceptance_rate_vector_fixed_eff[i] > sw_high_accep_rate) or (acceptance_rate_vector_fixed_eff_old[i] > sw_high_accep_rate and acceptance_rate_vector_fixed_eff[i] < sw_low_accep_rate):
                sw_adj_factor_fixed_eff[i] = 1 / np.sqrt(sw_adj_factor_fixed_eff[i])
            elif (acceptance_rate_vector_fixed_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_fixed_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_fixed_eff[i] < sw_low_accep_rate):
                sw_adj_factor_fixed_eff[i] = 2 / 3
            elif (acceptance_rate_vector_fixed_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_fixed_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_fixed_eff[i] > sw_high_accep_rate):
                sw_adj_factor_fixed_eff[i] = 3 / 2


    for i in range(len(sw_adj_factor_random_eff)):
        if equil_sweeps == delta_sw_adj_step:
            if (acceptance_rate_vector_random_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_random_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_random_eff[i] < sw_low_accep_rate) or (acceptance_rate_vector_random_eff[i] < sw_low_accep_rate):
                sw_adj_factor_random_eff[i] = 2 / 3
            elif (acceptance_rate_vector_random_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_random_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_random_eff[i] > sw_high_accep_rate) or (acceptance_rate_vector_random_eff[i] > sw_high_accep_rate):
                sw_adj_factor_random_eff[i] = 3 / 2
        else:
            if (acceptance_rate_vector_random_eff[i] <= sw_high_accep_rate and acceptance_rate_vector_random_eff[i] >= sw_low_accep_rate) or (step_width_random_eff[i] < sw_min_random_eff[i]) or (step_width_random_eff[i] > sw_max_random_eff[i]):
                sw_adj_factor_random_eff[i] = 1
            elif (acceptance_rate_vector_random_eff_old[i] < sw_low_accep_rate and acceptance_rate_vector_random_eff[i] > sw_high_accep_rate) or (acceptance_rate_vector_random_eff_old[i] > sw_high_accep_rate and acceptance_rate_vector_random_eff[i] < sw_low_accep_rate):
                sw_adj_factor_random_eff[i] = 1 / np.sqrt(sw_adj_factor_random_eff[i])
            elif (acceptance_rate_vector_random_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_random_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_random_eff[i] < sw_low_accep_rate):
                sw_adj_factor_random_eff[i] = 2 / 3
            elif (acceptance_rate_vector_random_eff_old[i] <= sw_high_accep_rate and acceptance_rate_vector_random_eff_old[i] >= sw_low_accep_rate and acceptance_rate_vector_random_eff[i] > sw_high_accep_rate):
                sw_adj_factor_random_eff[i] = 3 / 2

    for i in range(len(step_width_fixed_eff)):
        if step_width_fixed_eff[i] < sw_min_fixed_eff[i]:
            step_width_fixed_eff[i] = sw_min_fixed_eff[i]
        else:
            step_width_fixed_eff[i] = step_width_fixed_eff[i] * sw_adj_factor_fixed_eff[i]
    for i in range(len(step_width_random_eff)):
        if step_width_random_eff[i] < sw_min_random_eff[i]:
            step_width_random_eff[i] = sw_min_random_eff[i]
        else:
            step_width_random_eff[i] = step_width_random_eff[i] * sw_adj_factor_random_eff[i]
    acceptance_rate_vector_fixed_eff_old = acceptance_rate_vector_fixed_eff.copy()
    acceptance_rate_vector_random_eff_old = acceptance_rate_vector_random_eff.copy()
    acceptance_rate_vector_fixed_eff = np.zeros(len(acceptance_rate_vector_fixed_eff))
    acceptance_rate_vector_random_eff = np.zeros(len(acceptance_rate_vector_random_eff))
    acceptance_acum_vector_fixed_eff = np.zeros(len(acceptance_acum_vector_fixed_eff))
    acceptance_acum_vector_random_eff = np.zeros(len(acceptance_acum_vector_random_eff))

    return (step_width_fixed_eff , step_width_random_eff , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_rate_vector_fixed_eff_old , acceptance_rate_vector_random_eff_old , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff , sw_adj_factor_fixed_eff , sw_adj_factor_random_eff)



def mc_simulation(number_replica , initial_beta , initial_sigma_square , initial_rho , initial_gamma , data_values , y_values , variance_info , number_equilibration_sweeps , number_parameters , step_width_fixed_eff , step_width_random_eff , max_beta , min_beta , number_theta_vectors , number_decorrelation_sweeps , delta_sw_adj_step , sw_low_accep_rate , sw_high_accep_rate , name_simulation , seed_generator):
    np.random.seed(seed_generator)
    beta_vector_matrix = []
    gamma_vector_matrix = []
    sigma_square_vector_matrix = []
    file_output_beta = open(name_simulation + "/beta" , "w")
    file_output_gamma = open(name_simulation + "/gamma" , "w")
    file_output_sigma_square = open(name_simulation + "/sigma_square" , "w")
    for r in range(number_replica):
        #start of metropolis for each replica
        beta = initial_beta
        gamma = initial_gamma
        sigma_square = initial_sigma_square
        rho = initial_rho
        (V_matrix_inverse , log_determinant_V , dict_variance_order , dict_variance_block) = create_V_matrix(y_values , beta , data_values , variance_info , sigma_square)
        (hamiltonian , y_hat , y_diff) = hamiltonian_func(y_values , data_values , beta , V_matrix_inverse , log_determinant_V)
        variance_prediction = np.sum(np.square(y_diff)) / (len(y_diff) - 1)
        delta_variance_prediction = 0
        first_term = 0
        second_term = 0
        len_theta = len(beta) + len(sigma_square)


        #file_output_plot_test_ham = open("data_plot_test_ham_" + str(r) + "_" + name_simulation, "w")
        list_data_plot_test = []
        list_data_eq_beta = []
        list_data_eq_sigma = []
        list_data_eq_v_matrix_inverse = []
        list_data_eq_variance_prediction = []
        list_data_acceptance_betas = []
        list_data_acceptance_sigmas = []
        list_data_sw_betas = []
        list_data_sw_sigmas = []
        acceptance_rate_vector_fixed_eff = np.zeros(len(initial_beta))
        acceptance_rate_vector_random_eff = np.zeros(len(initial_sigma_square))
        acceptance_rate_vector_fixed_eff_old = acceptance_rate_vector_fixed_eff.copy()
        acceptance_rate_vector_random_eff_old = acceptance_rate_vector_random_eff.copy()
        acceptance_acum_vector_fixed_eff = np.zeros(len(initial_beta))
        acceptance_acum_vector_random_eff = np.zeros(len(initial_sigma_square))
        sweeps_number = 0
        acceptance_rate = 0
        sw_rate_fixed_eff = 1e-6
        sw_rate_random_eff = 1e-6
        sw_adj_factor_fixed_eff = np.ones(len(step_width_fixed_eff))
        sw_adj_factor_random_eff = np.ones(len(step_width_random_eff))
        step_width_fixed_eff_list = []
        step_width_random_eff_list = []
        for equil_sweeps in range(number_equilibration_sweeps):
            if equil_sweeps % delta_sw_adj_step == 0 and equil_sweeps != 0:
                list_data_acceptance_betas.append((acceptance_rate_vector_fixed_eff / acceptance_acum_vector_fixed_eff).reshape(len(acceptance_rate_vector_fixed_eff)))
                list_data_acceptance_sigmas.append((acceptance_rate_vector_random_eff / acceptance_acum_vector_random_eff).reshape(len(acceptance_rate_vector_random_eff)))
                list_data_sw_betas.append(step_width_fixed_eff.copy().reshape(len(step_width_fixed_eff)))
                list_data_sw_sigmas.append(step_width_random_eff.copy().reshape(len(step_width_random_eff)))
                (step_width_fixed_eff , step_width_random_eff , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_rate_vector_fixed_eff_old , acceptance_rate_vector_random_eff_old , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff , sw_adj_factor_fixed_eff , sw_adj_factor_random_eff) = step_width_adjustment(acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_rate_vector_fixed_eff_old , acceptance_rate_vector_random_eff_old , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff , sw_rate_fixed_eff , sw_rate_random_eff , step_width_fixed_eff , step_width_random_eff , sw_adj_factor_fixed_eff , sw_adj_factor_random_eff , sw_low_accep_rate , sw_high_accep_rate , delta_sw_adj_step , equil_sweeps , list_data_sw_betas , list_data_sw_sigmas)
            for theta_updaste in range(number_parameters):
                theta_selected = np.random.randint(0, len_theta)
                if theta_selected < len(beta):
                    delta_beta = (step_width_fixed_eff[theta_selected] * ((2 * np.random.uniform(0, 1)) - 1))
                    delta_rho = 0
                else:
                    delta_rho = (step_width_random_eff[theta_selected - len(beta)] * ((2 * np.random.uniform(0, 1)) - 1))
                    delta_beta = 0
                (beta , gamma , sigma_square , rho , hamiltonian , y_hat , V_matrix_inverse , log_determinant_V , y_diff , variance_prediction , delta_variance_prediction , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff) = mc_update(beta , gamma , sigma_square , rho , hamiltonian , y_values , y_hat , y_diff , variance_prediction , delta_variance_prediction , number_parameters , delta_beta , delta_rho , data_values , V_matrix_inverse , log_determinant_V , variance_info , max_beta , min_beta , theta_selected , dict_variance_order , dict_variance_block , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff)

            #----------------------------Plot stuff-------------------------------------------
            #file_output_plot_test_ham.write(str(hamiltonian) + "\t" + str(sweeps_number) + "\t" + str(first_term) + "\t" + str(second_term) + "\n")

            sweeps_number = sweeps_number + 1
            #------------------------------------------------------------------------------------
            list_data_plot_test.append([hamiltonian, sweeps_number, first_term, second_term])
            list_data_eq_beta.append(beta.copy().reshape(len(beta)))
            list_data_eq_sigma.append(sigma_square.copy().reshape(len(sigma_square)))
            list_data_eq_v_matrix_inverse.append(V_matrix_inverse.copy().reshape(len(V_matrix_inverse)))
            list_data_eq_variance_prediction.append(variance_prediction)

        np.savetxt(name_simulation + "/equilibration_beta" , np.array(list_data_eq_beta))
        np.savetxt(name_simulation + "/equilibration_sigma" , np.array(list_data_eq_sigma))
        np.savetxt(name_simulation + "/acceptance_fixed_eff" , np.array(list_data_acceptance_betas))
        np.savetxt(name_simulation + "/acceptance_random_eff" , np.array(list_data_acceptance_sigmas))
        np.savetxt(name_simulation + "/acceptance_sw_fixed_eff" , np.array(list_data_sw_betas))
        np.savetxt(name_simulation + "/acceptance_sw_random_eff" , np.array(list_data_sw_sigmas))
        np.savetxt(name_simulation + "/equilibration_hamiltonian" , np.array(list_data_plot_test))
        np.savetxt(name_simulation + "/equilibration_v_matrix_inverse" , np.array(list_data_eq_v_matrix_inverse))
        np.savetxt(name_simulation + "/equilibration_variance_prediction" , np.array(list_data_eq_variance_prediction))



        beta_vector_list = []
        gamma_vector_list = []
        sigma_square_vector_list = []
        for theta_acum in range(number_theta_vectors):
            for decorr_sweeps in range(number_decorrelation_sweeps):
                for theta_update in range(number_parameters):
                    theta_selected = np.random.randint(0, len_theta - 1)
                    if theta_selected < len(beta):
                        delta_beta = (step_width_fixed_eff[theta_selected] * ((2 * np.random.uniform(0, 1)) - 1))
                        delta_rho = 0
                    else:
                        delta_rho = (step_width_random_eff[theta_selected - len(beta)] * ((2 * np.random.uniform(0, 1)) - 1))
                        delta_beta = 0
                    (beta , gamma , sigma_square , rho , hamiltonian , y_hat , V_matrix_inverse , log_determinant_V , y_diff , variance_prediction , delta_variance_prediction , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff) = mc_update(beta , gamma , sigma_square , rho , hamiltonian , y_values , y_hat , y_diff , variance_prediction , delta_variance_prediction , number_parameters , delta_beta , delta_rho , data_values , V_matrix_inverse , log_determinant_V , variance_info , max_beta , min_beta , theta_selected , dict_variance_order , dict_variance_block , first_term , second_term , acceptance_rate_vector_fixed_eff , acceptance_rate_vector_random_eff , acceptance_acum_vector_fixed_eff , acceptance_acum_vector_random_eff)

                #-----------------Plot stuff--------------------------------------------------------------
                #file_output_plot_test_ham.write(str(hamiltonian) + "\t" + str(sweeps_number) + "\t" + str(first_term) + "\t" + str(second_term) + "\n")

                sweeps_number = sweeps_number + 1  # Plot stuff
                #-----------------------------------------------------------------------------------------
                list_data_plot_test.append([hamiltonian, sweeps_number, first_term, second_term])

            beta_vector_list.append(beta)
            gamma_vector_list.append(gamma)
            sigma_square_vector_list.append(sigma_square)

            for b_element in range(len(beta)):
                if b_element == len(beta) - 1:
                    file_output_beta.write(str(float(beta[b_element])))
                else:
                    file_output_beta.write(str(float(beta[b_element])) + " ")
            file_output_beta.write("|")
            for g_element in range(len(gamma)):
                if g_element == len(gamma) - 1:
                    file_output_gamma.write(str(float(gamma[g_element])))
                else:
                    file_output_gamma.write(str(float(gamma[g_element])) + " ")
            file_output_gamma.write("|")
            for sig_sq_element in range(len(sigma_square)):
                if sig_sq_element == len(sigma_square) - 1:
                    file_output_sigma_square.write(str(float(sigma_square[sig_sq_element])))
                else:
                    file_output_sigma_square.write(str(float(sigma_square[sig_sq_element])) + " ")
            file_output_sigma_square.write("|")



        beta_vector_matrix.append(beta_vector_list)
        gamma_vector_matrix.append(gamma_vector_list)
        sigma_square_vector_matrix.append(sigma_square_vector_list)
        file_output_beta.write("\n")
        file_output_gamma.write("\n")
        file_output_sigma_square.write("\n")
        #file_test.close()
        #file_output_plot_test_ham.close()
        np.savetxt(name_simulation + "/data_plot_test_ham_" + str(r) , np.array(list_data_plot_test))

    file_output_beta.close()
    file_output_gamma.close()
    file_output_sigma_square.close()
    #np.savetxt(name_simulation + "_beta" , np.array(beta_vector_matrix))
    #np.savetxt(name_simulation + "_gamma" , np.array(gamma_vector_matrix))
    #np.savetxt(name_simulation + "_sigma_square" , np.array(sigma_square_vector_matrix))
    print(acceptance_rate / ((number_equilibration_sweeps * number_parameters) + (number_theta_vectors * number_decorrelation_sweeps * number_parameters)))

    return (np.array(beta_vector_matrix) , np.array(gamma_vector_matrix))


def create_V_matrix(y_values , beta , data_values , variance_info , sigma_square):
    log_determinant = 0
    y_hat = np.dot(data_values, beta)
    y_diff = y_values - y_hat
    pred_diff_sum_square = np.sum(np.square(y_diff))
    variance_pred = pred_diff_sum_square / (len(y_diff) - 1)
    dict_variance_order = {}
    dict_variance_block = {}
    block_list = []
    for index , line in enumerate(open(variance_info)):
        data = line.split("|")
        index_processing = data[0].strip().split(" ")
        dict_variance_order[index] = list(map(lambda x: int(x), index_processing))
        for element in index_processing:
            dict_variance_block[int(element)] = index
        index_matrix_ini = int(index_processing[0]) #just one element, since all are the same
        X = data_values[index_matrix_ini]
        block_value = (np.matmul(X , X.T) * sigma_square[index_matrix_ini]) + variance_pred
        block_value_list = np.ones(len(index_processing)) * block_value
        block_list.extend(block_value_list)

    block_list_np = np.array(block_list)
    V_matrix_inverse = 1 / block_list_np
    V_matrix_log_determinant = np.sum(np.log(block_list))

    return (V_matrix_inverse , V_matrix_log_determinant , dict_variance_order , dict_variance_block)


def update_V_matrix(data_values , sigma_square , rho_sigma_change_index , V_matrix_inverse , variance_pred):
    index_matrix_ini = rho_sigma_change_index[0]
    X = data_values[index_matrix_ini]
    block_value = (np.matmul(X , X.T) * sigma_square[index_matrix_ini]) + variance_pred
    block_value_inverse = 1 / block_value
    for sigma_change_value in rho_sigma_change_index:
        V_matrix_inverse[sigma_change_value] = block_value_inverse

    block_list = 1 / V_matrix_inverse
    V_matrix_log_determinant = np.sum(np.log(block_list))

    return (V_matrix_inverse , V_matrix_log_determinant)











