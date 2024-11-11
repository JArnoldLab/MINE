import utilities as util
import numpy as np


list_sigma_selected = []
list_environmental_sigma_selected = []
list_pi_number = []
matrix_sigma = util.read_file("run_7/sigma_square")
matrix_sigma_average = util.average_theta_matrix(matrix_sigma)
matrix_sigma_average_individual = np.mean(matrix_sigma_average , axis = 0)
environmental_variance = util.get_initial_sigma_square("variance")

for index , line in enumerate(open("variance")):
    data = line.split("|")
    index_block = data[0].split(" ")
    index_selected = int(index_block[0])
    pi_number = data[1]
    if pi_number not in list_pi_number:
        list_sigma_selected.append(matrix_sigma_average_individual[index_selected])
        list_environmental_sigma_selected.append(environmental_variance[index_selected])
    list_pi_number.append(pi_number)

heritability = np.sum(list_sigma_selected) / (np.sum(list_sigma_selected) + np.sum(list_environmental_sigma_selected))#np.sum(matrix_sigma_average_individual) / (np.sum(matrix_sigma_average_individual) + np.sum(environmental_variance))

print(heritability)