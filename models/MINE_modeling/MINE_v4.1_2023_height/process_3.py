import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation = "run_7"




matrix_acceptance_betas = np.loadtxt(name_simulation + "/acceptance_fixed_eff")
matrix_acceptance_sigmas = np.loadtxt(name_simulation + "/acceptance_random_eff")
matrix_step_width_betas = np.loadtxt(name_simulation + "/acceptance_sw_fixed_eff")
matrix_step_width_sigmas = np.loadtxt(name_simulation + "/acceptance_sw_random_eff")





fig, ax = plt.subplots(2,2)
fig.tight_layout()
ax[0,0].plot(matrix_acceptance_betas)
ax.flat[0].set(ylabel = "Acceptance rate")
ax[0,0].set_title("Betas")



ax[0,1].plot(matrix_acceptance_sigmas)
ax[0,1].set_title("Sigmas")



ax[1,0].plot(matrix_step_width_betas)
ax.flat[2].set(ylabel = "Step width")
ax.flat[2].set(xlabel = "Sweeps x 1000")



ax[1,1].plot(matrix_step_width_sigmas)
ax.flat[3].set(xlabel = "Sweeps x 1000")




plt.savefig(name_simulation + "/step_width_combined" , bbox_inches = "tight")


