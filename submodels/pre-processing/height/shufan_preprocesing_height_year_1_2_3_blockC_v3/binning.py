import math
import numpy as np
import sys
import operator


input_file = input("Please enter vcf file name: ")
data_file = input("Please enter the phenotype file name: ")
min_number_bases = input("Please enter the minimum number of bases: ")
chromosome_number = input("Please enter the chromosome number: ")
total_snips = 0
number_accesion = 0
flag_accesion = 0
chromosome_dict = {'1':80880000 , '2':77740000 , '3':74390000 , '4':68660000 , '5':71850000 , '6':61280000 , '7':65510000 , '8':62690000 , '9':59420000 , '10':61230000}


genome_length = chromosome_dict[chromosome_number]

rec_snps_per_bin = genome_length / int(min_number_bases)

print("Chromosome size : " + str(genome_length))
print("Number of recommended snps per bin: " + str(rec_snps_per_bin))

number_snips  = input("Please enter the number of snps per bin: ")
number_snips_more = input("Please enter the number of snps to add in each bin if needed: ")
output_file = input("Please enter the name of the output file: ")


list_positions = []
for index, line in enumerate(open("../../../Sorghum_chromosomes/" + input_file)):
    if '/' in line:
        data = line.split("\t")
        total_snips = total_snips + 1
        if flag_accesion == 0:
            for i in range(len(data)):
                if '/' in data[i]:
                    number_accesion = number_accesion + 1
                    flag_accesion = 1

    if "POS" in line:
        labels = line.split("\t")
        for i in range(len(labels)):
            if ":" in labels[i]:
                name = labels[i].split(":")
                list_positions.append(name[0])



dict_accessions = {}
for index , line in enumerate(open(data_file)):
    if index != 0:
        data = line.split(",")
        if data[0] in dict_accessions:
            dict_accessions[data[0]].append(float(data[1].strip()))
        else:
            dict_accessions[data[0]] = []
            dict_accessions[data[0]].append(float(data[1].strip()))


dict_bins = {}
final_number_bases = 0



while final_number_bases == 0:

    number_bins = math.ceil(total_snips / int(number_snips))
    matrix_x = np.empty((number_accesion, int(number_bins)))


    count_bins = 0
    counter_row_data = 0
    bin_list = []
    final_number_bases = 1
    list_acum = []
    last_snip = 0
    first_snip = 0
    dict_bins_limits = {}

    for index , line in enumerate(open("../../../Sorghum_chromosomes/" + input_file)):
        if '/' in line:
            data = line.split("\t")
            list_line = []
            if counter_row_data % int(number_snips) == 0 and counter_row_data != 0:
                first_snip = int(data[1].strip())
                bases_length = first_snip - last_snip
                if bases_length < int(min_number_bases):
                    final_number_bases = 0
                    number_snips = int(number_snips) + int(number_snips_more)
                    break
                else:
                    dict_bins_limits[count_bins] = (last_snip , first_snip)
                    bin_list.append(bases_length)
                    matrix_x[:,count_bins] = list_acum
                    last_snip = first_snip
                    count_bins = count_bins + 1
                    list_acum = []

            for i in range(len(data)):
                if '/' in data[i]:
                    data_sum = data[i].strip().split("/")
                    if data_sum[0] != '.' and data_sum[1] != '.':
                        sum_1 = int(data_sum[0])
                        sum_2 = int(data_sum[1])
                        sum_total = sum_1 + sum_2
                        list_line.append(sum_total)
                    else:
                        list_line.append(0)
            if len(list_acum) == 0:
                list_acum = np.array(list_line)
            else:
                list_acum = list_acum + np.array(list_line)
            list_line = []
            counter_row_data = counter_row_data + 1

bin_list.append(genome_length - first_snip)
matrix_x[:,count_bins] = list_acum
dict_bins_limits[count_bins] = (first_snip , genome_length)


list_dry_weight_ordered = []
matrix_final = []
file_accession_labels = open(output_file + "_weight_labels" , "w")
file_accession_labels_variance = open(output_file + "_weight_variance" , "w")
index_matrix = 0
for i in range(len(list_positions)):
    key = list_positions[i]
    if key in dict_accessions:
        for values in dict_accessions[key]:
            matrix_final.append(matrix_x[i,:])
            list_dry_weight_ordered.append(values)
            file_accession_labels.write(key + "\n")
            file_accession_labels_variance.write(str(index_matrix) + " ")
            index_matrix = index_matrix + 1
        file_accession_labels_variance.write("|" + str(key) + "|" + str(np.var(dict_accessions[key])) + "\n")



matrix_final = np.array(matrix_final)
np.savetxt(output_file , matrix_final , fmt = '%1.2f')
np.savetxt(output_file + "_weight" , np.array(list_dry_weight_ordered))




print(matrix_final.shape)
print("Number of bins: " + str(len(bin_list)) + "\n")
for i in range(len(bin_list)):
    print("Bases in bin " + str(i) + ": " + str(bin_list[i]))

file_bin_location = open("chromosome_" + chromosome_number , "w")
for i in range(count_bins + 1):
    file_bin_location.write("bin_" + str(i + 1) + "," + str(dict_bins_limits[i][0]) + "," + str(dict_bins_limits[i][1]) + "\n")

file_bin_location.close()
