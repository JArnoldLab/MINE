chisq1 = readmatrix("MINE_v2_2023_2_disease_blockB/variance");
chisq2 = readmatrix("MINE_v2_2023_2_disease_blockC/variance");

variance = chisq1(:)
variance2 = chisq2(:)

figure(1)
h1 = histogram(variance)
h1.BinWidth = 0.001
xlabel('Chi-squared')
ylabel('Count')


figure(2)
h2 = histogram(variance2)
h2.BinWidth = 0.001
xlabel('Chi-squared')
ylabel('Count')


ham1 = readmatrix("MINE_v2_2023_2_disease_blockB/run_1/data_plot_test_ham_0");
ham2 = readmatrix("MINE_v2_2023_2_disease_blockC/run_1/data_plot_test_ham_0");

ham_hist = ham1(150000:300000,1)
ham_hist2 = ham2(150000:300000,1)

figure(3)
h3 = histogram(ham_hist)
h3.BinWidth = 0.01
xlabel('Chi-squared')
ylabel('Count')


figure(4)
h4 = histogram(ham_hist2)
h4.BinWidth = 0.01
xlabel('Chi-squared')
ylabel('Count')