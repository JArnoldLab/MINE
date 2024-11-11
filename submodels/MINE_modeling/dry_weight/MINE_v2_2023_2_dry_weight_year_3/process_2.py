import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation = "run_1"

matrix_hamiltonian = np.loadtxt(name_simulation + "/data_plot_test_ham_0")
matrix_betas_equilibrium = np.loadtxt(name_simulation + "/equilibration_beta")
matrix_acceptance_betas = np.loadtxt(name_simulation + "/acceptance")
matrix_step_width_betas = np.loadtxt(name_simulation + "/acceptance_sw")



plt.figure()
plt.plot(matrix_hamiltonian[:,0])
plt.xlabel("Sweeps")
plt.ylabel("Hamiltonian")
plt.title("Hamiltonian per sweeps")
plt.savefig(name_simulation + "/hamiltonian")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_hamiltonian[:,0]))
plt.xlabel("Sweeps")
plt.ylabel("Hamiltonian")
plt.title("Hamiltonian per sweeps")
plt.savefig(name_simulation + "/hamiltonian_log10")
plt.close("all")


plt.figure()
plt.plot(matrix_betas_equilibrium)
plt.xlabel("Sweeps")
plt.ylabel("Betas")
plt.title("Betas per sweeps")
plt.savefig(name_simulation + "/betas_eq")
plt.close("all")


plt.figure()
plt.plot(matrix_acceptance_betas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Acceptance rate")
plt.title("Betas acceptance rate per window of 1000")
plt.savefig(name_simulation + "/acceptance_betas")
plt.close("all")


plt.figure()
plt.plot(matrix_step_width_betas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Step width")
plt.title("Betas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_betas")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_step_width_betas))
plt.xlabel("Sweeps x 1000")
plt.ylabel("Log 10 Step width")
plt.title("Betas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_betas_log10")
plt.close("all")
