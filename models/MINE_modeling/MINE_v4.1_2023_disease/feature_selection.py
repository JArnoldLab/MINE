import numpy as np


def adjust_betas(betas_matrix , data_values , variance):
    variance_diag = np.diag(variance)
    x_t_mul_diag_inv = np.matmul(data_values.T , np.linalg.inv(variance_diag))
    a_matrix = np.matmul(x_t_mul_diag_inv , data_values)
    (eigval , eigvec) = np.linalg.eigh(a_matrix)
    rot = np.flip(eigvec , axis = 1)
    beta_star = np.matmul(rot.T , betas_matrix.T)
    for i in range(beta_star.shape[1]):
        for j in range(beta_star.shape[0]):
            if j > len(variance) - 1:
                beta_star[j , i] = 0

    beta_new = np.matmul(rot , beta_star)
    return(beta_new.T)