%compare blocks
clear all
phenotype = "height"
run = "run_1"
model = "MINE_v4.1_2023" % regular = MINE_v2_2023_2 , mixed = MINE_v4.1_2023
ini = 80000 % regular: 150000 mixed: disease= 15000 , dry_weight = 80000, height = 80000
fin = 1080000
model_directory_year_1 = strcat(phenotype , "/" , model , "_" , phenotype , "_year_1/" , run , "/" , "data_plot_test_ham_0")
model_directory_year_2 = strcat(phenotype , "/" , model , "_" , phenotype , "_year_2/" , run , "/" , "data_plot_test_ham_0")
model_directory_year_3 = strcat(phenotype , "/" , model , "_" , phenotype , "_year_3/" , run , "/" , "data_plot_test_ham_0")
chisq1 = readmatrix(model_directory_year_1);
chisq2 = readmatrix(model_directory_year_2);
chisq3 = readmatrix(model_directory_year_3);
figure(1)
histogram(chisq1(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Year 1 vs Year 2')
hold on
histogram(chisq2(ini:fin,1))
h1.normalization = 'probability'
figure(2)
histogram(chisq2(ini:fin,1))
h1.normalization = 'probability'
hold on
histogram(chisq3(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Year 2 vs Year 3')
% subplot
figure(3)
subplot(3,1,1) 
histogram(chisq1(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Year 1')
xlim([-750 , -400])
h1.normalization = 'probability'
subplot(3,1,2)
histogram(chisq2(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Year 2')
xlim([-750 , -400])
h1.normalization = 'probability'
subplot(3,1,3)
histogram(chisq3(ini:fin,1))
xlabel('Hamiltonian')
ylabel('Count')
title('Year 3')
xlim([-750 , -400])
h1.normalization = 'probability'
figure(4)
histogram(chisq1(ini:fin,1))
h1.normalization = 'probability'
xlabel('Hamiltonian')
ylabel('Count')
title('Year 1 vs Year 3')
hold on
histogram(chisq3(ini:fin,1))
h1.normalization = 'probability'