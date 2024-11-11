%compare blocks
clear all
phenotype = "disease"
run = "run_1"
model = "MINE_v4.1_2023"% regular = MINE_v2_2023_2 , mixed = MINE_v4.1_2023
ini = 15001 % regular: 150000 mixed: disease= 15000 , dry_weight = 120000, height = 80000
fin = 1015000
model_directory_year_2 = strcat(phenotype , "/" , model , "_" , phenotype , "_year_2/" , run , "/" , "data_plot_test_ham_0")
model_directory_year_3 = strcat(phenotype , "/" , model , "_" , phenotype , "_year_3/" , run , "/" , "data_plot_test_ham_0")
chisq2 = readmatrix(model_directory_year_2);
chisq3 = readmatrix(model_directory_year_3);
figure(1)
histogram(chisq2(ini:fin,1))
h1.normalization = 'probability'
hold on
histogram(chisq3(ini:fin,1))
h1.normalization = 'probability'
xlabel('Chi-squared')
ylabel('Count')
title('Year 2 vs Year 3')
% subplot
figure(2)
subplot(3,1,1)
histogram(chisq2(ini:fin,1))
xlabel('Chi-squared')
ylabel('Count')
title('Year 2')
h1.normalization = 'probability'
subplot(3,1,2)
histogram(chisq3(ini:fin,1))
xlabel('Chi-squared')
ylabel('Count')
title('Year 3')
h1.normalization = 'probability'