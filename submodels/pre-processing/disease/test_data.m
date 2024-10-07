


data_year_2 = readmatrix("shufan_preprocesing_disease_year_2/disease_log10.csv");
data_year_3 = readmatrix("shufan_preprocesing_disease_year_3/disease_log10.csv");
figure(1)
h1 = histogram(data_year_2(:,2))
%h1.Normalization = 'probability'
xlabel('disease score')
ylabel('Probability')
title('year 2 vs year 3')
hold on
h2 = histogram(data_year_3(:,2))
%h2.Normalization = 'probability'


figure(2)
subplot(3,1,1) 
h3 = histogram(data_year_2(:,2))
xlabel('disease variance')
ylabel('Probability')
title('year 2')
%h3.Normalization = 'probability'
subplot(3,1,2)
h4 = histogram(data_year_3(:,2))
xlabel('disease variance')
ylabel('Probability')
title('year 3')
%h4.Normalization = 'probability'



var_year_2 = readmatrix("shufan_preprocesing_disease_year_2/chr1_weight_variance");
var_year_3 = readmatrix("shufan_preprocesing_disease_year_3/chr1_weight_variance");
figure(3)
h5 = histogram(var_year_2)
%h5.Normalization = 'probability'
h5.BinWidth = 0.001
xlabel('disease variance')
ylabel('Probability')
title('variance year 2 vs variance year 3')
hold on
h6 = histogram(var_year_3)
%h6.Normalization = 'probability'
h6.BinWidth = 0.001


figure(4)
subplot(3,1,1) 
h7 = histogram(var_year_2)
xlabel('disease variance')
ylabel('Probability')
title('variance year 2')
%h7.Normalization = 'probability'
h7.BinWidth = 0.001
subplot(3,1,2)
h8 = histogram(var_year_3)
xlabel('disease variance')
ylabel('Probability')
title('variance year 3')
%h8.Normalization = 'probability'
h8.BinWidth = 0.001
