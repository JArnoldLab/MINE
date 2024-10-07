


dict_gene_id = {}

for index , line in enumerate(open("summary_genes_height.csv")):
    data = line.strip().split("\t")
    dict_gene_id[data[0]] = data[2:]

file_output = open("table_genes_height_temp.csv" , "w")
for index , line in enumerate(open("table_genes_height.csv")):
    data = line.strip().split("\t")
    if data[0] in dict_gene_id:
        if len(data) > 1:
            file_output.write(line.strip() + "\t" + dict_gene_id[data[0]][0] + "\t" + dict_gene_id[data[0]][1] + "\t" + dict_gene_id[data[0]][2] + "\n")
        else:
            file_output.write(line.strip() + "\t\t" + dict_gene_id[data[0]][0] + "\t" + dict_gene_id[data[0]][1] + "\t" + dict_gene_id[data[0]][2] + "\n")

    else:
        file_output.write(line)

file_output.close()