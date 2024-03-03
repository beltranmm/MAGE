
%% simulate gene expression profile

% load MAGE functions
f = MAGE_functions;

% set profile parameters
numGene = 1000;
numSamplesPerCondition = 4;
mu = 7;
SD = 1;
DErate = 0.1;

% generate expression profile
[profile,DEind] = f.simData(numGene,numSamplesPerCondition,mu,SD,DErate);

% display mean expression
figure;
scatter(mean(profile(:,1:4),2),mean(profile(:,5:8),2))
xlabel('Condition A (expression level)')
ylabel('Condition B (expression level)')
title('Original mean expression')


%% filter out low expression

% remove genes with expr < 1 in < 6 samples
[profile,fltrind] = f.fltr(profile,1,6);

numGene = size(profile,1);

% update DE gene indices
DEind = DEind(fltrind);

% display filtered mean expression
figure;
scatter(mean(profile(:,1:4),2),mean(profile(:,5:8),2))
xlabel('Condition A (expression level)')
ylabel('Condition B (expression level)')
title('Filtered mean expression')


%% find aberrantly expressed genes using MAGE

[~,OS] = MAGE(profile(:,1:4),profile(:,5:8),100,5,1,0.9);

%% find differentially expressed genes based on fold-change

FC = mean(profile(:,1:4),2) - mean(profile(:,5:8),2);

%% find differentially expressed genes based on 2-sample t-test

pVal = zeros(numGene,1);
for i = 1 : numGene
    [~,pVal(i)] = ttest2(profile(i,1:4),profile(i,5:8));
end

%% Analysis: ROC/AUC

% ROC
[TPR_mage,FPR_mage] = f.ROC(OS,DEind);
[TPR_fc,FPR_fc] = f.ROC(abs(FC),DEind);
[TPR_p,FPR_p] = f.ROC(1-pVal,DEind);

figure;
plot(FPR_mage,TPR_mage)
hold on
plot(FPR_fc,TPR_fc)
plot(FPR_p,TPR_p)
plot([0,1],[0,1],'--k')
hold off
legend({'MAGE','FC','t-test','random'},'Location','southeast')
xlabel('false-positive rate')
ylabel('true-positive rate')
title('ROC analysis')

% AUC
disp(' ')
disp('AUC')
disp('---------')
disp("MAGE: " + string(round(f.AUC(FPR_mage,TPR_mage),2)))
disp("FC: " + string(round(f.AUC(FPR_fc,TPR_fc),2)))
disp("p-value: " + string(round(f.AUC(FPR_p,TPR_p),2)))
disp("random: 0.50")




