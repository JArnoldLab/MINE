from PyPDF2 import PdfReader
import os


gene_file = "genes_Kresovich.csv"
pdf_directory = "pdfs"
list_genes = []

for index , line in enumerate(open(gene_file)):
    if line.strip() != "":
        list_genes.append(line.strip())

(path , dirs , files) = next(os.walk(pdf_directory))
file_output = open("summary_" + gene_file , "w")

for gene in list_genes:
    for file in files:
        flag = 0
        reader = PdfReader(pdf_directory + "/" + file)
        for i in range(len(reader.pages)):
            if gene.upper() in reader.pages[i].extract_text().upper():
                flag = 1
                break
        if flag == 1:
            file_output.write(gene + "\t" + file + "\n")

file_output.close()
