from biomart import BiomartServer


def genes_retrieval_sorghum(region_begin , region_end , chromosome):
    server = BiomartServer("http://phytozome-next.jgi.doe.gov/biomart")

    server.verbose = True

    phytozome = server.datasets['phytozome']

    response = phytozome.search({
    'filters': {'organism_id': '454' , 'region_end' : [region_end] , 'region_name' : [chromosome] , 'region_start' : [region_begin]},
    'attributes': ['gene_name1', 'gene_description' , 'chr_name1']
    })

    return response.text



def search_location_bin(index_bin , phenotype):
    index_total = 0
    flag_found = 0
    i = 1
    while flag_found == 0 and i <= 10:
        for index , line in enumerate(open("chromosomes/chromosome_" + str(i))):
            if index_total == index_bin:
                data = line.split(",")
                begin = int(data[1])
                end = int(data[2].strip())
                if i < 10:
                    chromosome = "chr0" + str(i)
                else:
                    chromosome = "chr" + str(i)
                flag_found = 1
                break
            index_total = index_total + 1
        i = i + 1
    return(begin , end , chromosome)


def find_genes(bin_indexes , phenotype , name_simulation):
    file_output = open("genes_" + phenotype + "_" + name_simulation , "w")
    file_output_err = open("errors_" + phenotype + "_" + name_simulation , "w")
    file_output.write("Gene_name\tDescription\tChromosome\n")
    for index in bin_indexes:
        (begin , end , chromosome) = search_location_bin(index , phenotype)
        try:
            genes = genes_retrieval_sorghum(begin , end , chromosome)
            file_output.write("--------------\tBin " + str(index + 1) + "\t------------------\n")
            file_output.write(genes)
        except:
            file_output_err.write("Error in bin " + str(index + 1) + "\n")
    file_output.close()
    file_output_err.close()

