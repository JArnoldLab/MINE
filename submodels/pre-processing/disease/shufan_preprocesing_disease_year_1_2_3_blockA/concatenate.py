import numpy as np

file_names = input("Please enter the files to concatenate. eg. chr1,chr2,chr3 : ")
output_file = input("Please enter the name of the final matrix: ")


files = file_names.split(",")
acum_matrix = []

for file_chromosome in files:
    matrix_temp = []
    for index , line in enumerate(open(file_chromosome)):
        data = line.strip().split(" ")
        data = list(map(lambda x:float(x) , data))
        matrix_temp.append(data)

    if len(acum_matrix) == 0:
        acum_matrix = np.array(matrix_temp)
    else:
        acum_matrix = np.concatenate((acum_matrix , np.array(matrix_temp)) , axis = 1)


norms_matrix = np.linalg.norm(acum_matrix , axis = 0)
acum_matrix = acum_matrix / norms_matrix
print(acum_matrix.shape)
np.savetxt(output_file , acum_matrix , fmt = '%1.4f')
