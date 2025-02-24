Running the programs:

MINE

#Pre-processing
For pre-processing run the script binning.py for each chromosome to obtain the design matrix X from the VCF file on Sorghum. The use is the following:
python binning.py
Once each gene matrix is generated, use concatenate.py to form the final design matrix X. The use is the following:
python concatenate.py
#Modeling and MINE procedure
Fill the CONFIG file with the options you want, either for modeling or the MINE procedure. The use is the following:
python main_file.py CONFIG
To obtain the figures proving the modeling went well, use the script process_2.py, edit the file providing the name of the run you wrote previously in the config file. The use is the following:
python process_2.py
The folders with prefix MINE_v2 represent the linear model, and the folders with prefix MINE_v4.1 represent the mixed linear model.
#Gene finding
Copy the chromosomal region limits from pre-processing in a folder called chromosomes. Use the main.py file, edit it according to the name of the run and other files you are providing.  This will retrieve the genes in the chromosomal bins that survived feature selection in the ensemble run. The use is the following:
python main.py
#Pdf scanning on genes
Put the paper pdfs you want your genes scanned on in the folder pdfs. Use the script main_extraction.py, edit the file according to the files you have. The use is the following:
python main_extraction.py
#Additional files
There are a few script files such as those to obtain figures in matlab from partitioned data in years or blocks. They are in the submodels folder.



Description of files

Genes_pdf_extraction

This file contains all the genes extracted from the literature and a comparison of those in our study for: (i) dry weight; (ii) height; (iii) disease as well as those extracted from Kresovich data for dry weight. The python program is main_extraction.py

MINE_modeling

These files containing the ensemble modeling and MINE calculations for: (i) dry weight; (ii) height; (iii) disease as well as for Kresovich data on dry weight.  The analysis includes the software for generating the design matrix X and heritability, feature selection, and the main program for the ensemble run. Each  MINE file contains up to 11 MCMC runs under varied conditions for year and block. Varied plots of statistics like stepwidth versus sweep are included. Versions with the number of 4.1 include the option of the mixed linear model; Versions with the label v2 include the option of the linear model. The main program for doing the ensemble modeling is main_file.py.  There are also Matlab files for displaying the results of an ensemble run or MINE calculation in some subfolders.

Data vectors that have names like disease, height, and dry weight.  The X matrix has a name like matrix_dry_weight_log20.

Within each folder for a run under MINE modeling there is one file that executes the ensemble method and MINE.  It is called main_file.py.  That is the main program to execute the Markov Chain Monte Carlo(MCMC) Experiment . For example, we have under:

MINE/models/MINE_modeling/MINE_v4.1_2023_dry_weight/

The output of this run is in folder run_14.  This run involves a mixed linear model (v4.1).



Gene_finder

These files contain the results of running the Gene_extraction program with Biomart for: (i) dry weight; (ii) height; and (iii) disease.  The program for doing the gene finding is called main.py.

Pre-processing 

These folders contain the generation of the X-matrix and sample variances from the VCF files for: (i) dry weight; (ii) height; (iii) disease. The program binning.py carries out pre-processing for individual chromosomes.

Submodels

This folder contains substitute models for: (i) dry weight; (ii) height; (iii) disease for the main model involving no block or year effects involving all data from years 1-3 and Kresovich data.  For example the substitutes includes models that allow year effects or block effects. There are also Matlab files to display the results.

Take the following folder as an example:

MINE_V4.1_2023_dry_weight_year_2_no_blockA

This run sets aside block A and runs an MCMC experiments with all other year 2 blocks on dry weight.  The main program is main_file.py.  The output is in run_1 folder.

VCF files

These VCF files for sorghum are too large for GitHub.  They can be found on a dropbox folder with address given in the manuscript in the data and code availability section.
The programs can be run so long as there is a process design matrix X derived from the VCF files.



