import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation_1 = "run_1"
name_simulation_2 = "run_1"
name_simulation_3 = "run_1"
name_model_1 = "MINE_v4.1_2023_dry_weight_year_3"
name_model_2 = "MINE_v4.1_2023_height_year_3"
name_model_3 = "MINE_v4.1_2023_disease_year_3"


matrix_hamiltonian_1 = np.loadtxt("dry_weight/" + name_model_1 + "/" + name_simulation_1 + "/data_plot_test_ham_0")
matrix_hamiltonian_2 = np.loadtxt("height/" + name_model_2 + "/" + name_simulation_2 + "/data_plot_test_ham_0")
matrix_hamiltonian_3 = np.loadtxt("disease/" + name_model_3 + "/" + name_simulation_3 + "/data_plot_test_ham_0")



fig , (ax1,ax2,ax3) = plt.subplots(1,3)
fig.tight_layout()
ax1.plot(matrix_hamiltonian_1[:,0])
ax1.set(ylabel = "Hamiltonian")
ax1.set_title("Dry weight")
ax2.plot(matrix_hamiltonian_2[:,0])
ax2.set(xlabel = "Sweeps")
ax2.set_title("Height")
ax3.plot(matrix_hamiltonian_3[:,0])
ax3.set_title("Disease")
plt.savefig("hamiltonian_mixed_year_3" , bbox_inches = "tight")




