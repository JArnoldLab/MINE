%compare blocks
clear all
chisqyr1 = readmatrix('data_plot_test_ham_0_run_1_year1');
chisqyr2 = readmatrix('data_plot_test_ham_0_run_1_year2');
figure(1)
histogram(chisqyr1(500000:1040000,1))
xlabel('Chi-squared')
ylabel('Probability')
h1.normalization = 'probability'
hold on
histogram(chisqyr2(500000:1040000,1))
h1.normalization = 'probability'
xlabel('Chi-squared')
ylabel('Probability')
% subplot
figure(3)
subplot(2,1,1) 
histogram(chisqyr1(500000:1040000,1))
xlabel('chi-square')
ylabel('probability')
h1.normalization = 'probability'
subplot(2,1,2)
histogram(chisqyr2(500000:1040000,1))
h1.normalization = 'probability'
xlabel('chi-square')
ylabel('probability')