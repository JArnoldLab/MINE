import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation_1 = "run_1"
name_simulation_2 = "run_1"
name_simulation_3 = "run_1"
name_model_1 = "MINE_v2_2023_2_dry_weight_year_2"
name_model_2 = "MINE_v2_2023_2_height_year_2"
name_model_3 = "MINE_v2_2023_2_disease_year_2"



 
matrix_betas_equilibrium_1 = np.loadtxt("dry_weight/" + name_model_1 + "/" + name_simulation_1 + "/equilibration_beta")
matrix_betas_equilibrium_2 = np.loadtxt("height/" + name_model_2 + "/" + name_simulation_2 + "/equilibration_beta")
matrix_betas_equilibrium_3 = np.loadtxt("disease/" + name_model_3 + "/" + name_simulation_3 + "/equilibration_beta")




fig , (ax1,ax2,ax3) = plt.subplots(1,3)
ax1.plot(matrix_betas_equilibrium_1)
ax1.set(ylabel = "Betas")
ax1.set_title("Dry weight")
ax2.plot(matrix_betas_equilibrium_2)
ax2.set(xlabel = "Sweeps")
ax2.set_title("Height")
ax3.plot(matrix_betas_equilibrium_3)
plt.savefig("betas_year_2")


