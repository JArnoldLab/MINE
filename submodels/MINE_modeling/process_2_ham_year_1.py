import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation_1 = "run_1"
name_simulation_2 = "run_1"
name_simulation_3 = "run_1"
name_model_1 = "MINE_v2_2023_2_dry_weight_year_1"
name_model_2 = "MINE_v2_2023_2_height_year_1"
name_model_3 = "MINE_v2_2023_2_disease_year_1"

matrix_hamiltonian_1 = np.loadtxt("dry_weight/" + name_model_1 + "/" + name_simulation_1 + "/data_plot_test_ham_0")
matrix_hamiltonian_2 = np.loadtxt("height/" + name_model_2 + "/" + name_simulation_2 + "/data_plot_test_ham_0")
#matrix_hamiltonian_3 = np.loadtxt("disease/" + name_model_3 + "/" + name_simulation_3 + "/data_plot_test_ham_0")
matrix_betas_equilibrium_1 = np.loadtxt("dry_weight/" + name_model_1 + "/" + name_simulation_1 + "/equilibration_beta")
matrix_betas_equilibrium_2 = np.loadtxt("height/" + name_model_2 + "/" + name_simulation_2 + "/equilibration_beta")
#matrix_betas_equilibrium_3 = np.loadtxt("disease/" + name_model_3 + "/" + name_simulation_3 + "/equilibration_beta")





fig , (ax1,ax2) = plt.subplots(1,2)
ax1.plot(np.log10(matrix_hamiltonian_1[:,0]))
ax1.set(ylabel = "Log10 Hamiltonian")
ax1.set(xlabel = "Sweeps")
ax1.set_title("Dry weight")
ax2.plot(np.log10(matrix_hamiltonian_2[:,0]))
ax2.set(xlabel = "Sweeps")
ax2.set_title("Height")
#ax3.plot(np.log10(matrix_hamiltonian_3[:,0]))
#ax3.set_title("Disease")
plt.savefig("hamiltonian_year_1")


fig , (ax1,ax2) = plt.subplots(1,2)
ax1.plot(matrix_betas_equilibrium_1)
ax1.set(ylabel = "Betas")
ax1.set(xlabel = "Sweeps")
ax1.set_title("Dry weight")
ax2.plot(matrix_betas_equilibrium_2)
ax2.set(xlabel = "Sweeps")
ax2.set_title("Height")
#ax3.plot(np.log10(matrix_hamiltonian_3[:,0]))
#ax3.set_title("Disease")
plt.savefig("betas_year_1")

