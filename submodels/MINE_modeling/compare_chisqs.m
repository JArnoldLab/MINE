%compare blocks
clear all
phenotype = "height"
run = "run_1"
model = "MINE_v4.1_2023" % regular = MINE_v2_2023_2 , mixed = MINE_v4.1_2023
ini = 80000 % regular: 150000 , disease = 300000 mixed: disease= 15000 , dry_weight = 80000, height = 80000
fin = 1080000
model_directory_blockA = strcat(phenotype , "/" , model , "_" , phenotype , "_blockA/" , run , "/" , "data_plot_test_ham_0")
model_directory_blockB = strcat(phenotype , "/" , model , "_" , phenotype , "_blockB/" , run , "/" , "data_plot_test_ham_0")
model_directory_blockC = strcat(phenotype , "/" , model , "_" , phenotype , "_blockC/" , run , "/" , "data_plot_test_ham_0")
chisqA = readmatrix(model_directory_blockA);
chisqB = readmatrix(model_directory_blockB);
chisqC = readmatrix(model_directory_blockC);
figure(1)
histogram(chisqA(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Block A vs. Block B')
hold on
histogram(chisqB(ini:fin,1))
h1.normalization = 'probability'
figure(2)
histogram(chisqB(ini:fin,1))
h1.normalization = 'probability'
hold on
histogram(chisqC(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Block B vs Block C')
% subplot
figure(3)
subplot(3,1,1) 
histogram(chisqA(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Block A')
h1.normalization = 'probability'
subplot(3,1,2)
histogram(chisqB(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Block B')
h1.normalization = 'probability'
subplot(3,1,3)
histogram(chisqC(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Block C')
h1.normalization = 'probability'
figure(4)
histogram(chisqA(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Block A vs. Block C')
hold on
histogram(chisqC(ini:fin,1))
h1.normalization = 'probability'