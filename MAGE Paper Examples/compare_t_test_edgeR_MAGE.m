% compare t-test and edgeR results
%% setup

% select data directory
dataDir = strcat(pwd,'\workspaces\edgeR\');

% load contour functions
f = DEG_contour_functions;

%% edgeR vs t-test (gt3)

[gene_R,FC_R,pVal_R] = import_edgeR(strcat(dataDir,'benchmark_edgeR_gt3.csv'));
[gene,FC_t,FDR_t,OS,FDR_M] = import_tt_MAGE(strcat(dataDir,'MAGE and t-test results for breast cancer gT3 treatment profile.csv'));

% filter and sort for common genes
[ind,ind_R] = geneFilter(gene,gene_R);
gene = gene(ind);
FC_t = FC_t(ind);
FDR_t = FDR_t(ind);
OS = OS(ind);
FDR_M = FDR_M(ind);
gene_R = gene_R(ind_R);
FC_R = FC_R(ind_R);
pVal_R = pVal_R(ind_R);

ind = geneSort(gene,gene_R);
gene_R = gene_R(ind);
FC_R = FC_R(ind);
pVal_R = pVal_R(ind);



% load profile data
[profile,geneName,sampleName] = f.load_gt3();
sampleClass = ["control";"treatment"];


% t-test p-values
pVal_t = zeros(numel(geneName),1);
for i = 1 : numel(pVal_t)
    [~,pVal_t(i)] = ttest2(profile(i,strcmp(sampleName,sampleClass(1))),...
        profile(i,strcmp(sampleName,sampleClass(2))));
end


% filter and sort for common genes
[ind,ind_R] = geneFilter(geneName,gene_R);
geneName = geneName(ind);
pVal_t = pVal_t(ind);
gene_R = gene_R(ind_R);
pVal_R = pVal_R(ind_R);

ind = geneSort(geneName,gene_R);
gene_R = gene_R(ind);
pVal_R = pVal_R(ind);

figure;
scatter(pVal_t,pVal_R)
xlabel('t-test: (p-value)')
ylabel('edgeR: (p-value)')

figure;
scatter(-log10(pVal_t),-log10(pVal_R),'filled')
xlabel('t-test: (-log_{10}p-value)')
ylabel('edgeR: (-log_{10}p-value)')
if max(xlim) > max(ylim)
    ylim([0,max(xlim)])
else
    xlim([0,max(ylim)])
end

%clear profile geneName sampleName sampleClass i

%% edgeR vs t-test (mTOR)

[gene_R,FC_R,pVal_R] = import_edgeR(strcat(dataDir,'benchmark_edgeR_mtor.csv'));
[gene,FC_t,FDR_t,OS,FDR_M] = import_tt_MAGE(strcat(dataDir,'MAGE and t-test results for mTOR KO mouse profile.csv'));

% filter and sort for common genes
[ind,ind_R] = geneFilter(gene,gene_R);
gene = gene(ind);
FC_t = FC_t(ind);
FDR_t = FDR_t(ind);
OS = OS(ind);
FDR_M = FDR_M(ind);
gene_R = gene_R(ind_R);
FC_R = FC_R(ind_R);
pVal_R = pVal_R(ind_R);

ind = geneSort(gene,gene_R);
gene_R = gene_R(ind);
FC_R = FC_R(ind);
pVal_R = pVal_R(ind);



% load mTOR data
[profile,geneName,sampleName] = f.load_mTOR();
sampleClass = ["control";"treatment"];


% t-test p-values
pVal_t = zeros(numel(geneName),1);
for i = 1 : numel(pVal_t)
    [~,pVal_t(i)] = ttest2(profile(i,strcmp(sampleName,sampleClass(1))),...
        profile(i,strcmp(sampleName,sampleClass(2))));
end


% filter and sort for common genes
[ind,ind_R] = geneFilter(geneName,gene_R);
geneName = geneName(ind);
pVal_t = pVal_t(ind);
gene_R = gene_R(ind_R);
pVal_R = pVal_R(ind_R);

ind = geneSort(geneName,gene_R);
gene_R = gene_R(ind);
pVal_R = pVal_R(ind);

figure;
scatter(pVal_t,pVal_R)
xlabel('t-test: (p-value)')
ylabel('edgeR: (p-value)')

figure;
scatter(-log10(pVal_t),-log10(pVal_R),'filled')
xlabel('t-test: (-log_{10}p-value)')
ylabel('edgeR: (-log_{10}p-value)')
if max(xlim) > max(ylim)
    ylim([0,max(xlim)])
else
    xlim([0,max(ylim)])
end

%clear profile geneName sampleName sampleClass i
    

%% functions
% import t-test and MAGE results
function [geneName_t,FC_t,FDR_t,OS,FDR_M] = import_tt_MAGE(fileName)

    benchmark = readmatrix(fileName,'OutputType','string');

    geneName_t = benchmark(:,1);
    geneName_M = benchmark(:,1);
    OS = str2double(benchmark(:,2));
    FDR_M = str2double(benchmark(:,3));
    FC_t = str2double(benchmark(:,4));
    FDR_t = str2double(benchmark(:,5));
    clear benchmark
end


% import edgeR results
function [geneName_edgeR,FC,pVal] = import_edgeR(fileName)
    
    benchmark_edgeR = readmatrix(fileName,'OutputType','string');

    geneName_edgeR = benchmark_edgeR(:,1);
    FC = str2double(benchmark_edgeR(:,2));
    CPM = str2double(benchmark_edgeR(:,3));
    Lratio = str2double(benchmark_edgeR(:,4));
    pVal = str2double(benchmark_edgeR(:,5));
    F = str2double(benchmark_edgeR(:,6));
    FpVal = str2double(benchmark_edgeR(:,7));
    clear benchmark_edgeR
end

% filter by gene names
function [indA,indB] = geneFilter(listA,listB)
    indA = ismember(listA,listB);
    indB = ismember(listB,listA);
end

% sort by gene names
function ind = geneSort(listA,listB)
    [~,ind] = ismember(listA,listB);
end


