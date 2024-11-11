chisq1 = readmatrix("variance");

variance = chisq1(:,3)

figure(1)
h1 = histogram(variance)
xlabel('Chi-squared')
ylabel('Count')