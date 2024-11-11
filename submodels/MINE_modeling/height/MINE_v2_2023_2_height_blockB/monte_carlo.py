import numpy as np
import matplotlib.pyplot as plt
import random
import time
import utilities as util

plt.style.use('ggplot')



def hamiltonian_func(beta , gamma , data_values , y_values , variance):
    prediction = np.dot(data_values , beta).reshape(len(data_values) , 1)
    fraction = variance * (np.square(y_values - prediction))
    sum_experiments = np.sum(fraction)

    return (sum_experiments , prediction)

def hamiltonian_func_fast(delta_beta , data_values , y_values , y_previous , variance_rec , beta_selected):
    delta_y = (data_values[:,beta_selected] * delta_beta).reshape(y_values.shape)
    y_new = y_previous + delta_y
    hamiltonian = np.sum(np.square(y_values - y_new) * variance_rec)
    return(hamiltonian , y_new)




def ensemble(delta_hamiltonian):
    return np.exp(-delta_hamiltonian / 2)


def mc_update(beta , gamma , hamiltonian , y_previous , number_betas , step_width , data_values , y_values , variance , max_beta , min_beta , acceptance_rate_vector , acceptance_acum_vector):
    acceptance = 0
    beta_new = beta.copy()
    gamma_new = gamma.copy()
    len_theta = len(beta)# + len(gamma)
    theta_selected = random.randint(0, len_theta - 1)
    beta_new[theta_selected] = beta_new[theta_selected] + (step_width[theta_selected] * ((2 * np.random.uniform(0, 1)) - 1))  # in the general case I would alter on the chosen variable
    acceptance_acum_vector[theta_selected] = acceptance_acum_vector[theta_selected] + 1
    delta_beta = beta_new[theta_selected] - beta[theta_selected]
    (hamiltonian_new , y_new) = hamiltonian_func_fast(delta_beta , data_values , y_values , y_previous , variance , theta_selected)

    delta_hamiltonian = hamiltonian_new - hamiltonian

    if delta_hamiltonian < 0:
        ratio = 2.0
    else:
        ratio = ensemble(delta_hamiltonian)

    if ratio >= 1:
        acceptance = 1
        acceptance_rate_vector[theta_selected] = acceptance_rate_vector[theta_selected] + 1
    else:
        v = np.random.uniform(0, 1)
        if v < ratio:
            acceptance = 1
            acceptance_rate_vector[theta_selected] = acceptance_rate_vector[theta_selected] + 1


    if acceptance == 1:
        return (beta_new , gamma_new , hamiltonian_new , y_new , acceptance_rate_vector , acceptance_acum_vector)
    else:
        return (beta , gamma , hamiltonian , y_previous , acceptance_rate_vector , acceptance_acum_vector)


def step_width_adjustment(acceptance_rate_vector , acceptance_rate_vector_old , acceptance_acum_vector , sw_rate , step_width , sw_adj_factor , sw_low_accep_rate , sw_high_accep_rate , delta_sw_adj_step , equil_sweeps):
    acceptance_rate_vector = acceptance_rate_vector / acceptance_acum_vector
    sw_min = sw_rate * np.mean(step_width)
    sw_max = (1 / sw_rate) * np.mean(step_width)
    for i in range(len(sw_adj_factor)):
        if equil_sweeps == delta_sw_adj_step:
            if (acceptance_rate_vector_old[i] <= sw_high_accep_rate and acceptance_rate_vector_old[i] >= sw_low_accep_rate and acceptance_rate_vector[i] < sw_low_accep_rate) or (acceptance_rate_vector[i] < sw_low_accep_rate):
                sw_adj_factor[i] = 2 / 3
            elif (acceptance_rate_vector_old[i] <= sw_high_accep_rate and acceptance_rate_vector_old[i] >= sw_low_accep_rate and acceptance_rate_vector[i] > sw_high_accep_rate) or (acceptance_rate_vector[i] > sw_high_accep_rate):
                sw_adj_factor[i] = 3 / 2
        else:
            if (acceptance_rate_vector[i] <= sw_high_accep_rate and acceptance_rate_vector[i] >= sw_low_accep_rate) or (step_width[i] < sw_min) or (step_width[i] > sw_max):
                sw_adj_factor[i] = 1
            elif (acceptance_rate_vector_old[i] < sw_low_accep_rate and acceptance_rate_vector[i] > sw_high_accep_rate) or (acceptance_rate_vector_old[i] > sw_high_accep_rate and acceptance_rate_vector[i] < sw_low_accep_rate):
                sw_adj_factor[i] = 1 / np.sqrt(sw_adj_factor[i])
            elif (acceptance_rate_vector_old[i] <= sw_high_accep_rate and acceptance_rate_vector_old[i] >= sw_low_accep_rate and acceptance_rate_vector[i] < sw_low_accep_rate):
                sw_adj_factor[i] = 2 / 3
            elif (acceptance_rate_vector_old[i] <= sw_high_accep_rate and acceptance_rate_vector_old[i] >= sw_low_accep_rate and acceptance_rate_vector[i] > sw_high_accep_rate):
                sw_adj_factor[i] = 3 / 2

    step_width = step_width * sw_adj_factor
    acceptance_rate_vector_old = acceptance_rate_vector.copy()
    acceptance_rate_vector = np.zeros(len(acceptance_rate_vector))
    acceptance_acum_vector = np.zeros(len(acceptance_acum_vector))

    return (step_width , acceptance_rate_vector , acceptance_rate_vector_old , acceptance_acum_vector , sw_adj_factor)




def mc_simulation(number_replica , initial_beta , initial_gamma , data_values , y_values , variance , number_equilibration_sweeps , number_betas , step_width , max_beta , min_beta , number_theta_vectors , number_decorrelation_sweeps , delta_sw_adj_step , name_simulation):
    np.random.seed(int(time.time()))
    beta_vector_matrix = []
    gamma_vector_matrix = []
    acceptance_rate_vector = np.zeros(len(initial_beta))
    acceptance_rate_vector_old = acceptance_rate_vector.copy()
    acceptance_acum_vector = np.zeros(len(initial_beta))
    sweeps_number = 0
    sw_low_accep_rate = 0.30
    sw_high_accep_rate = 0.70
    sw_rate = 1e-6
    sw_adj_factor = np.ones(len(step_width))
    file_output_beta = open(name_simulation + "/beta" , "w")
    file_output_gamma = open(name_simulation + "/gamma" , "w")
    for r in range(number_replica):
        #start of metropolis for each replica
        beta = initial_beta
        gamma = initial_gamma
        (hamiltonian , y_previous) = hamiltonian_func(beta , gamma, data_values , y_values , variance)


        #file_output_plot_test_ham = open("data_plot_test_ham_" + str(r) + "_" + name_simulation, "w")
        list_data_hamiltonian = []
        list_data_eq_beta = []
        list_data_acceptance_betas = []
        list_data_sw_betas = []
        for equil_sweeps in range(number_equilibration_sweeps):
            if equil_sweeps % delta_sw_adj_step == 0 and equil_sweeps != 0:
                list_data_acceptance_betas.append((acceptance_rate_vector / acceptance_acum_vector).reshape(len(acceptance_rate_vector)))
                list_data_sw_betas.append(step_width.reshape(len(step_width)))
                (step_width , acceptance_rate_vector , acceptance_rate_vector_old , acceptance_acum_vector , sw_adj_factor) = step_width_adjustment(acceptance_rate_vector , acceptance_rate_vector_old , acceptance_acum_vector , sw_rate , step_width , sw_adj_factor , sw_low_accep_rate , sw_high_accep_rate , delta_sw_adj_step , equil_sweeps)
            for theta_updaste in range(number_betas):
                (beta , gamma , hamiltonian , y_previous , acceptance_rate_vector , acceptance_acum_vector) = mc_update(beta , gamma , hamiltonian , y_previous , number_betas , step_width , data_values , y_values , variance , max_beta , min_beta , acceptance_rate_vector , acceptance_acum_vector)

            #----------------------------Plot stuff-------------------------------------------
            #file_output_plot_test_ham.write(str(hamiltonian) + "\t" + str(sweeps_number) + "\n")

            sweeps_number = sweeps_number + 1
            #------------------------------------------------------------------------------------
            list_data_hamiltonian.append([hamiltonian , sweeps_number])
            list_data_eq_beta.append(beta.reshape(len(beta)))

        np.savetxt(name_simulation + "/equilibration_beta" , np.array(list_data_eq_beta))
        np.savetxt(name_simulation + "/acceptance" , np.array(list_data_acceptance_betas))
        np.savetxt(name_simulation + "/acceptance_sw" , np.array(list_data_sw_betas))
        np.savetxt(name_simulation + "/equilibration_hamiltonian" , np.array(list_data_hamiltonian))
        list_data_eq_beta = []
        list_data_acceptance_betas = []
        list_data_sw_betas = []

        beta_vector_list = []
        gamma_vector_list = []
        for theta_acum in range(number_theta_vectors):
            for decorr_sweeps in range(number_decorrelation_sweeps):
                for theta_update in range(number_betas):
                    (beta , gamma , hamiltonian , y_previous , acceptance_rate_vector , acceptance_acum_vector) = mc_update(beta , gamma , hamiltonian , y_previous , number_betas , step_width , data_values , y_values , variance , max_beta , min_beta , acceptance_rate_vector , acceptance_acum_vector)

                #----------------------------Plot stuff-------------------------------------------------
                #file_output_plot_test_ham.write(str(hamiltonian) + "\t" + str(sweeps_number) + "\n")

                sweeps_number = sweeps_number + 1
                #-----------------------------------------------------------------------------------------
                list_data_hamiltonian.append([hamiltonian , sweeps_number])

            beta_vector_list.append(beta)
            gamma_vector_list.append(gamma)
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


        beta_vector_matrix.append(beta_vector_list)
        gamma_vector_matrix.append(gamma_vector_list)
        file_output_beta.write("\n")
        file_output_gamma.write("\n")
        #file_test.close()
        #file_output_plot_test_ham.close()
        np.savetxt(name_simulation + "/data_plot_test_ham_" + str(r) , np.array(list_data_hamiltonian))

    file_output_beta.close()
    file_output_gamma.close()
    print(acceptance_rate_vector / ((number_equilibration_sweeps * number_betas) + (number_theta_vectors * number_decorrelation_sweeps * number_betas)))

    return (np.array(beta_vector_matrix) , np.array(gamma_vector_matrix))



