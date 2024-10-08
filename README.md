# MINE

#Pre-processing


For pre-processing run the script binning.py for each chromosome. The use is the following:

python binning.py

Once each gene matrix is generated, use concatenate.py to form the final design matrix. The use is the following:

python concatenate.py


#Modeling and MINE procedure

Fill the CONFIG file with the options you want, either for modeling or the MINE procedure. The use is the following:

python main_file.py CONFIG

To obtain the figures proving the modeling went well, use the script process_2.py, edit the file providing the name of the run you wrote previously in the config file. The use is the following:

python process_2.py

The folders with prefix MINE_v2 represent the linear model, and the folders with prefix MINE_v4.1 represent the mixed linear model.


#Gene finding

Copy the chromosomal region limits from pre-processing in a folder called chromosomes. Use the main.py file, edit it according to the name of the run and other files you are providing. The use is the following:

python main.py


#Pdf scanning on genes

Put the paper pdfs you want your genes scanned on in the folder pdfs. Use the script main_extraction.py, edit the file according to the files you have. The use is the following:

python main_extraction.py


#Additional files

There are a few script files such as those to obtain figures in matlab from partitioned data in years or blocks. They are in the submodels folder.



























