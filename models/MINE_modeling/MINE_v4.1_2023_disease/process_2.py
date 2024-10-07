import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

name_simulation = "run_1"

matrix_hamiltonian = np.loadtxt(name_simulation + "/data_plot_test_ham_0")
matrix_betas_equilibrium = np.loadtxt(name_simulation + "/equilibration_beta")
matrix_sigmas_equilibrium = np.log10(np.loadtxt(name_simulation + "/equilibration_sigma"))
matrix_acceptance_betas = np.loadtxt(name_simulation + "/acceptance_fixed_eff")
matrix_acceptance_sigmas = np.loadtxt(name_simulation + "/acceptance_random_eff")
matrix_step_width_betas = np.loadtxt(name_simulation + "/acceptance_sw_fixed_eff")
matrix_step_width_sigmas = np.loadtxt(name_simulation + "/acceptance_sw_random_eff")
matrix_v_matrix_inverse = np.loadtxt(name_simulation + "/equilibration_v_matrix_inverse")
matrix_environmental_sigma = np.loadtxt(name_simulation + "/equilibration_variance_prediction")



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
plt.ylabel("Log 10 Hamiltonian")
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
plt.plot(matrix_sigmas_equilibrium)
plt.xlabel("Sweeps")
plt.ylabel("Log10 sigmas")
plt.title("Sigmas per sweeps")
plt.savefig(name_simulation + "/sigmas_eq")
plt.close("all")


plt.figure()
plt.plot(matrix_acceptance_betas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Acceptance rate")
plt.title("Betas acceptance rate per window of 1000")
plt.savefig(name_simulation + "/acceptance_betas")
plt.close("all")


plt.figure()
plt.plot(matrix_acceptance_sigmas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Acceptance rate")
plt.title("Sigmas acceptance rate per window of 1000")
plt.savefig(name_simulation + "/acceptance_sigmas")
plt.close("all")


plt.figure()
plt.plot(matrix_step_width_betas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Step width")
plt.title("Betas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_betas")
plt.close("all")


plt.figure()
plt.plot(matrix_step_width_sigmas)
plt.xlabel("Sweeps x 1000")
plt.ylabel("Step width")
plt.title("Sigmas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_sigmas")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_step_width_betas))
plt.xlabel("Sweeps x 1000")
plt.ylabel("Log 10 Step width")
plt.title("Betas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_betas_log10")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_step_width_sigmas))
plt.xlabel("Sweeps x 1000")
plt.ylabel("Log 10 Step width")
plt.title("Sigmas step width per window of 1000")
plt.savefig(name_simulation + "/step_width_sigmas_log10")
plt.close("all")


plt.figure()
plt.plot(matrix_v_matrix_inverse)
plt.xlabel("Sweeps")
plt.ylabel("V matrix inverse")
plt.title("V matrix inverse per sweeps")
plt.savefig(name_simulation + "/v_matrix_inverse")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_v_matrix_inverse))
plt.xlabel("Sweeps")
plt.ylabel("Log 10 V matrix inverse")
plt.title("V matrix inverse per sweeps")
plt.savefig(name_simulation + "/v_matrix_inverse_log10")
plt.close("all")


plt.figure()
plt.plot(1 / matrix_v_matrix_inverse)
plt.xlabel("Sweeps")
plt.ylabel("V matrix")
plt.title("V matrix per sweeps")
plt.savefig(name_simulation + "/v_matrix")
plt.close("all")


plt.figure()
plt.plot(np.log10(1 / matrix_v_matrix_inverse))
plt.xlabel("Sweeps")
plt.ylabel("Log 10 V matrix")
plt.title("V matrix per sweeps")
plt.savefig(name_simulation + "/v_matrix_log10")
plt.close("all")


plt.figure()
plt.plot(matrix_environmental_sigma)
plt.xlabel("Sweeps")
plt.ylabel("Environmental sigma")
plt.title("Environmental sigma per sweeps")
plt.savefig(name_simulation + "/variance_prediction")
plt.close("all")


plt.figure()
plt.plot(np.log10(matrix_environmental_sigma))
plt.xlabel("Sweeps")
plt.ylabel("Log 10 Environmental sigma")
plt.title("Environmental sigma per sweeps")
plt.savefig(name_simulation + "/variance_prediction_log10")
plt.close("all")
