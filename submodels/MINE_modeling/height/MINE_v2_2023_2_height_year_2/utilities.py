import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import operator
import os
from datetime import date

plt.style.use('ggplot')





def read_data_plot(file_name):
    theta_matrix = []
    list_sweeps = []
    for index , line in enumerate(open(file_name)):
        data = line.split("\t")
        theta_vector = data[0].split(",")
        theta_vector = list(map(lambda x: float(x.strip()) , theta_vector))
        theta_matrix.append(theta_vector)
        list_sweeps.append(int(data[2].strip()))
    return (np.array(theta_matrix) , list_sweeps)

def contour_plot(data_values , recommended_data_d):
    x = []
    y = []
    z = []
    log_det = []
    triangles_alt = []
    for i in range(len(recommended_data_d)):
        (x_comp , y_comp , z_comp) = data_values[recommended_data_d[i][1]][0]
        log_det_comp = recommended_data_d[i][0]
        x.append(x_comp)
        y.append(y_comp)
        z.append(z_comp)
        log_det.append(log_det_comp)
    triangles = mtri.Triangulation(x , y)
    plt.tricontourf(triangles , log_det)
    plt.xlabel("u1 coordinate")
    plt.ylabel("u2 coordinate")
    plt.title("Contour plot based on u vector coordinates and volume value")
    plt.show()


def plot_chi_sq_sweeps(file ,  name_simulation):
    steps = []
    chi_sq_e = []
    chi_sq_d = []
    for index , line in enumerate(open(file)):
        if index != 0:
            data = line.split("|")
            chi_sq_d.append(math.log10(float(data[0])))
            chi_sq_e.append(math.log10(float(data[1])))
            steps.append(int(data[2].strip()))

    plt.figure()
    plt.plot(steps , chi_sq_d)
    plt.xlabel("Steps")
    plt.ylabel("Log(det(D-matrix))")
    plt.title("Hamiltonian trajectory - Covariance criteria")
    plt.savefig("chi-squared_d" + name_simulation , bbox_inches = 'tight')
    plt.show()

    plt.figure()
    plt.plot(steps, chi_sq_e)
    plt.xlabel("Steps")
    plt.ylabel("Log(det(E-matrix))")
    plt.title("Hamiltonian trajectory - Correlation criteria")
    plt.savefig("chi-squared_e" + name_simulation, bbox_inches='tight')
    plt.show()



def plot_chi_sq_sweeps_mine(file ,  name_simulation):
    steps = []
    chi_sq_e = []
    chi_sq_d = []
    for index , line in enumerate(open(file)):
        if index != 0 and index <= 2000:
            data = line.split("|")
            chi_sq_d.append(float(data[0]))
            chi_sq_e.append(float(data[1]))
            steps.append(int(data[2].strip()))

    plt.figure()
    plt.plot(steps , chi_sq_d)
    plt.xlabel("Steps")
    plt.ylabel("Log(det(D-matrix))")
    plt.title("Hamiltonian trajectory - Covariance criteria")
    plt.savefig("chi-squared_d_" + name_simulation , bbox_inches = 'tight')
    plt.show()

    plt.figure()
    plt.plot(steps, chi_sq_e)
    plt.xlabel("Steps")
    plt.ylabel("Log(det(E-matrix))")
    plt.title("Hamiltonian trajectory - Correlation criteria")
    plt.savefig("chi-squared_e_" + name_simulation, bbox_inches='tight')
    plt.show()


def hist_chi_squared(ham_result , title):
    plt.hist(ham_result)
    plt.title("MC sample Hamiltonian" + " " + title)
    plt.show()

def hist_chi_squared_difference(ham_result_1 , ham_result_2):
    difference = ham_result_1 - ham_result_2
    plt.hist(abs(difference))
    plt.title("MC sample Hamiltonian difference (beta and beta star)")
    plt.show()


def plot_betas_sweeps(file , number_equilibration_sweeps , name_simulation):
    betas = []
    sweeps = []
    for index , line in enumerate(open(file)):
        data = line.split("\t")
        data_betas = data[0].split(" ")
        data_betas = list(map(lambda x: float(x[1 : len(x) - 1]) , data_betas))
        betas.append(data_betas)
        sweeps.append(int(data[1].strip()))

    plt.figure()
    plt.plot(sweeps , betas)
    plt.xlabel("Equilibration sweeps")
    plt.ylabel("Beta parameters")
    plt.title("Beta parameters at every equilibration sweep")
    plt.savefig("betas_sweeps_" + name_simulation , bbox_inches = 'tight')
    plt.show()




def plot_log_chi_sq_sweeps(chi_sq , sweeps , number_equilibration_sweeps , name_simulation):
    plt.figure()
    chi_sq = np.log10(chi_sq)
    plt.plot(sweeps , chi_sq)
    plt.axvline(x = number_equilibration_sweeps , color = 'blue')

    plt.xlabel("Eql sweeps: 1 to " + str(number_equilibration_sweeps) + " , Acc: " + str(number_equilibration_sweeps) + " to " + str(len(sweeps)))
    plt.ylabel("Log chi-Squared")
    plt.title("Trajectory Log(Chi-Squared) in Eq and Acc stages")
    plt.savefig("log_chi-squared_" + name_simulation , bbox_inches = 'tight')
    plt.show()


def plot_theta_sweeps(theta_matrix , sweeps , number_equilibration_sweeps , name_simulation , real_theta):
    plt.figure()
    plt.plot(sweeps , theta_matrix[:,0] , label = "theta_1")
    plt.plot(sweeps , theta_matrix[:,1] , label = "theta_2")
    plt.plot(sweeps , theta_matrix[:,2] , label = "theta_3")
    plt.axvline(x = number_equilibration_sweeps , color = 'blue')
    for theta in real_theta:
        plt.axhline(y = theta , color = 'orange')

    plt.xlabel("Eql sweeps: 1 to " + str(number_equilibration_sweeps) + " , Acc: " + str(number_equilibration_sweeps) + " to " + str(len(sweeps)))
    plt.ylabel("Thetas")
    plt.title("Trajectory Thetas in Eq and Acc stages")
    plt.legend()
    plt.savefig("thetas_" + name_simulation , bbox_inches = 'tight')
    plt.show()

def plot_hist_thetas(theta_matrix , name_simulation , max_theta , number_bins):
    plt.figure()
    plt.hist(theta_matrix[:,0] , range = (15 , 18) , histtype = 'barstacked' , rwidth = 0.85 , label = 'Theta_1')
    plt.xlabel("Theta 1 values")
    plt.ylabel("count (theta_1)")
    plt.title("Histogram: Theta 1 in Eql and Acc stages")
    plt.savefig("hist_theta_1_" + name_simulation , bbox_inches = 'tight')
    plt.show()

    plt.figure()
    plt.hist(theta_matrix[:, 1] , range = (50 , 52) , histtype = 'barstacked', rwidth = 0.85, label = 'Theta_2')
    plt.xlabel("Theta 2 values")
    plt.ylabel("count (theta_2)")
    plt.title("Histogram: Theta 2 in Eql and Acc stages")
    plt.savefig("hist_theta_2_" + name_simulation , bbox_inches = 'tight')
    plt.show()

    plt.figure()
    plt.hist(theta_matrix[:, 2] , range = (80 , 82) , histtype = 'barstacked', rwidth = 0.85,label = 'Theta_3')
    plt.xlabel("Theta 3 values")
    plt.ylabel("count (theta_3)")
    plt.title("Histogram: Theta 3 in Eql and Acc stages")
    plt.savefig("hist_theta_3_" + name_simulation , bbox_inches = 'tight')
    plt.show()


def plot_hist_chi_sq(chi_sq , name_simulation , max_theta , number_bins):
    plt.figure()
    plt.hist(chi_sq , range = (-10 , 10) , histtype = 'barstacked', rwidth = 0.85, label = 'Chi-squared')
    plt.xlabel("Chi-Squared values")
    plt.ylabel("count(chi-squared)")
    plt.title("Histogram: Chi-Squared in Eql and Acc stages")
    plt.savefig("hist_chi_sq_" + name_simulation , bbox_inches = 'tight')
    plt.show()


def plot_log_det_rank(file , name_simulation , number_points_plot):
    list_data_d_matrix = []
    list_data_e_matrix = []
    for index, line in enumerate(open(file)):
        if index != 0:
            data = line.split("|")
            list_data_d_matrix.append((data[0], float(data[1])))
            list_data_e_matrix.append((data[0], float(data[2].strip())))
    list_data_d_matrix = sorted(list_data_d_matrix, key=operator.itemgetter(1), reverse=True)
    list_data_e_matrix = sorted(list_data_e_matrix, key=operator.itemgetter(1), reverse=True)
    (_ , log_det_e) = zip(*list_data_e_matrix)
    (_ , log_det_d) = zip(*list_data_d_matrix)

    plt.figure()
    plt.plot(list(range(number_points_plot)) , log_det_e[:number_points_plot])
    plt.xlabel("Rank")
    plt.ylabel("Log(det(E))")
    plt.title("Top " + str(number_points_plot) + " tuples correlation score")
    plt.savefig("rank_e_" + name_simulation, bbox_inches='tight')
    plt.show()

    plt.figure()
    plt.plot(list(range(number_points_plot)), log_det_d[:number_points_plot])
    plt.xlabel("Rank")
    plt.ylabel("Log(det(D))")
    plt.title("Top " + str(number_points_plot) + " tuples covariance score")
    plt.savefig("rank_d_" + name_simulation, bbox_inches='tight')
    plt.show()

def plot_box_plot_each_run(file_prefix_boxplot , number_initial_simulation , number_simulations , mine_option):
    list_data_box_plot = []
    list_positions = []
    for i in range(number_simulations):
        for index , line in enumerate(open("run_" + str(number_initial_simulation + i) + file_prefix_boxplot)):
            data = line.split("|")
            if mine_option == 1:
                log_det = float(data[2].strip())
            else:
                log_det = float(data[1].strip())
            list_data_box_plot.append(log_det)
    plt.figure()
    plt.boxplot(list_data_box_plot , positions = [1])
    plt.xlabel("Set of experiments")
    if mine_option == 1:
        plt.ylabel("Log(Det(E))")
    else:
        plt.ylabel("Log(Det(D))")
    plt.title("Log Det in " + str(number_simulations) + " expriments")
    plt.savefig("log_det_n_experiments_box_plot_" + str(mine_option) , bbox_inches = 'tight')
    plt.show()

    plt.figure()
    for value in list_data_box_plot:
        plt.scatter(1 , value)
    plt.xlabel("Set of experiments")
    if mine_option == 1:
        plt.ylabel("Log(Det(E))")
    else:
        plt.ylabel("Log(Det(D))")
    plt.title("Log Det in " + str(number_simulations)+ " experiments")
    plt.savefig("log_det_n_experiments_" + str(mine_option)  , bbox_inches = 'tight')
    plt.show()


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

def read_config_file(config_file):
    dict_config_file = {}
    for index , line in enumerate(open(config_file)):
        if len(line.split("=")) == 2:
            data = line.split("=")
            dict_config_file[data[0].strip().upper()] = data[1].strip()
    path = dict_config_file['NAME_SIMULATION']
    os.system("rm -rf " + path)
    os.mkdir(path)
    os.system("cp " + config_file + " " + path)

    return dict_config_file


def average_theta_matrix(theta_matrix):
    (number_replica , number_theta_vectors_acum , number_theta_parameters) = theta_matrix.shape
    matrix_avg = np.zeros((number_theta_vectors_acum , number_theta_parameters))
    for i in range(number_replica):
        matrix_avg = matrix_avg + theta_matrix[i]
    matrix_avg = matrix_avg / number_replica
    return matrix_avg



def do_plots(name_simulation , number_replica , number_steps_annealing , number_equilibration_sweeps):

    #file_rank = name_simulation + "_labels_results"
    #file_rank_sq = name_simulation + "_chi_sq"
    file_prefix_boxplot = "_analysis_labels_results"
    #file_chi_sq = "data_plot_test_ham_" + str(number_replica) + "_" + name_simulation
    #file_betas = "data_plot_test_beta_" + str(number_replica) + "_" + name_simulation
    #plot_log_det_rank(file_rank , name_simulation , number_points_plot)
    #plot_chi_sq_sweeps(file_chi_sq , name_simulation)
    #plot_chi_sq_sweeps_mine(file_rank_sq , name_simulation)
    #plot_betas_sweeps(file_betas , number_equilibration_sweeps , name_simulation)
    plot_box_plot_each_run(file_prefix_boxplot , 26 , 10 , 0)


