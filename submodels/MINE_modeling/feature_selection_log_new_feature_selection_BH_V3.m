clear all
% This version includes an extra plot where the Y and Yhat are normalized 
% for the difference variances of the accessions
% start = sum of the number of bins from all earlier chromosomes + 1
% last  = start + number of bins of current chromosome.
start = 1048
last  = 1350
% get the ensemble
%uses projection method
%uses bayesian confidence interval for feature selection
%uses Benjamini-Hochberg for feature selection
betas = readmatrix('beta_star_dry_weight');
X = readmatrix('matrix_dry_weight_log20');
Y = readmatrix('dry_weight');
sig2 = readmatrix('variance')
% keep only the beta columns (markers) that are significant
n = 1000
p = 2748
% sort betas
betaranked = sort(betas,1);
% find betas where 2.5% and 97.5% don't include 0
for j=1:p
    nokeep(j) = 0;
end
for j=1:p
    if betaranked(25,j) < 0 & betaranked(975,j) > 0
        nokeep(j) = 1;
    end
end
% create the keepers
keep = ones(1,p) - nokeep;
dfr = sum(keep);
% do Benjamini-Hochberg
% calculate mean and standard error for each beta
alpha = 0.05
meanbeta = mean(betas);
sdbeta = std(betas);
%sebeta = sdbeta/1000^.5;
sebeta = sdbeta;
z = meanbeta./sebeta;
pd = makedist('normal',0,1);
zabs = abs(z);
pr = cdf(pd,zabs,'upper');
pr = 2*pr;
[sortp,key] = sort(pr);
for i=1:p
    cutoff(i) = i;
end
cutoff = alpha*cutoff/p;
% determine what to keep
for i=1:p
    keep2(i) = 0;
    if sortp(i) < cutoff(i);
        keep2(i) = 1;
    end
end
figure(1)
plot(1:i,sortp,1:i,cutoff)
xlabel('rank of P-value')
ylabel('contribution to P-value')
title('Benjamini-Hochberg Feature Selection')
% use key to retrieve the position of keep2
for i=1:p
    keep3(i) = 0;
end
    for i=1:p
        if keep2(i) > 0
        keep3(key(i)) = 1;
        end
    end
    % Benjamini-Hochberg is above
    % Now combine Bayesian Interval(keep) with BH(keep3).
product = keep3.*keep;
%product = keep  % Bayesian interval only
%product = keep3  % Benjamini-Hochberg only
j = 0
for i=1:p
    if product(i) == 1
        j=j+1;
        betakeep3(j) = betas(i);
    end
end
nkeep3 = j
%
% create X matrix that only keeps the keepers
jj = 0
for j=1:p
      %if product(j) == 1
 %    if keep(j) == 1 & j >= 1048 & j <= 1350 % choose a chromosome
      if product(j) == 1 & j >= start & j <= last % choose a chromosome
        jj = jj + 1;
        Xkeep(:,jj) = X(:,j);
        betakeep(:,jj) = betas(:,j);
        betakeeploc(:,jj) = j;
        pkeep(jj) = pr(j)
    end
end
nkeep= jj
for j=1:p
    chrom5(j) = 0;
end
for j=start:last
   chrom5(j) = 1;
end
meanbeta = mean(betakeep,1);
meanbeta = meanbeta';
yhat = Xkeep*meanbeta;
ehat = Y - yhat;
RSS = yhat'*yhat;
TSS = Y'*Y;
ESS = ehat'*ehat;
figure(2)
normplot(ehat)
xlabel('residual')
ylabel('ranked residual from normal')
figure(3)
scatter(yhat,ehat)
xlabel('predicted Y')
ylabel('residual')
R2 = RSS/(RSS+ESS)  % uses meanbeta
%pull out chromosome 5
jj = 0
betakeep5loc = zeros(1,dfr);
for j=1:nkeep
    if betakeeploc(j) >= start & betakeeploc(j) <= last
        jj = jj + 1;
        betakeep5loc(jj) = betakeeploc(j);
    end
end
% check Monte Carlo
figure(4)
chisq = readmatrix('data_plot_test_ham_0_run_1');
plot(chisq(:,2),chisq(:,1))
ylim([0 10000])
xlabel('sweeps')
ylabel('chisq')
figure(5)
scatter(key,-log(sortp))
xlabel('chromosome bin')
ylabel('-logP')
title('-logP for all chromosome bins')
figure(6)
scatter(betakeeploc,-log(pkeep))
title('chromosome regions feature selected')
xlabel('chromosome bin')
ylabel('-logP')
% do Manhattan plot with all filters applied using product
j = 0
for i=1:p
    %if product(i) == 1
    if product(i) == 1 & i >= start & i <= last % choose a chromosome
        j = j + 1
        key3(j) = i
        manh(j) = -log(pkeep(j))
    end
end
figure(7)
scatter(key3,manh)
% make correction for unequal variances in residuals
sig = sig2.^.5
yhatnorm = yhat./sig
Ynorm = Y./sig
ehatnorm = Ynorm - yhatnorm
figure(8)
scatter(yhatnorm,ehatnorm)
title('Y and Yhat normalized')
xlabel('predicted Y')
ylabel('residual')
RSSnorm = yhatnorm'*yhatnorm
ESSnorm = ehatnorm'*ehatnorm
TSSnorm = Ynorm'*Ynorm