%% setup

% select data directory
dataDir = strcat(pwd,'\workspaces\');

% load contour functions
f = DEG_contour_functions;


profile_complete = zeros();

%% load profile: breast cancer mouse mTOR KO treatment/control
load(strcat(dataDir,'breast_cancer_mouse_mTOR_delta_NCBI.mat'));

geneName = geneName_WT100;
clear geneName_WT100 geneName_WT101 geneName_WT105 geneName_KO107 geneName_KO111 geneName_KO115
sampleName = ["control";"control";"control";"treatment";"treatment";"treatment"];
profile = [profile_WT100 profile_WT101 profile_WT105 profile_KO107 profile_KO111 profile_KO115];
clear profile_WT100 profile_WT101 profile_WT105 profile_KO107 profile_KO111 profile_KO115

% remove all zero FPKM
ind = find(sum(profile,2) > 0);
profile = profile(ind,:);
geneName = geneName(ind);


profile = f.logNorm(profile,1,0);

% create legend, group, and colors
lgnd = unique(sampleName);

grp = zeros(numel(sampleName),1);
for i = 1 : numel(sampleName)
    grp(i) = find(strcmp(sampleName(i),lgnd));
end
clrset = hot(numel(lgnd));
grpclr = clrset;

title1 = "breast cancer (mouse mTOR KO)";

%[final_profile, filteredInd] = fltr(profile,level,sampleMin)
%f.filterAnalysis(profile,[0.1,1,3,5,8],[1,2,3,6]);
[~, ind] = f.fltr(profile,0.1,5);
profile = profile(ind,:);
geneName = geneName(ind);

profile_complete = profile;
geneName_complete = geneName;

%% simple 1D/2D (treatment/control)

%MAGE parameters
gridDensity = 100;
numContours = 5;
targetContainment = 0.95;
removeHighLowExpr = 0;

% add noise
gaussianNoise = false; % not implemented yet
increaseSTDX = 0;
increaseSTDY = 0;

% run topo
runTopo = false;

% run MAGE
runMAGE = true;

% use edgeR results
edgeR = true;

% Classification Thresholds
FCcutoff = 0.5;
pValcutoff = 0.05;
outlierScoreCutoff = 0.2;
outlierDistScoreCutoff = -10;

% pathway enrichment
exportForPE = 1;
useTopNgenes = 200;

%trackGene = ["EGR1";"ATF3";"DDIT3";"PTPRO";"GDF15";"MAFF";"WIPI1"];
trackGene = [];

sampleClass = ["control";"treatment"];

% analysis
HC = 0;


% sample permutation FDR
permFDR = true;
modAvoidFDROverestimate = false;

% save output to csv
saveToTable = true;
tableName = 'breast_cancer_gt3_MAGE_analysis.csv';


%--------------------------------------------------------------------------
profile = profile_complete;
geneName = geneName_complete;

% increase noise
if increaseSTDX > 0
    ind = find(strcmp(sampleName,sampleClass(1)));
    for j = 1 : numel(ind)
        for i = 1 : numel(geneName)
            profile(i,ind(j)) = profile(i,ind(j)) + increaseSTDX*randn(1,1);
        end
    end
end
if increaseSTDY > 0
    ind = find(strcmp(sampleName,sampleClass(2)));
    for j = 1 : numel(ind)
        for i = 1 : numel(geneName)
            profile(i,ind(j)) = profile(i,ind(j)) + increaseSTDY*randn(1,1);
        end
    end
end

meanExpr_T = mean(profile(:,strcmp(sampleName,sampleClass(1))),2);
meanExpr_C = mean(profile(:,strcmp(sampleName,sampleClass(2))),2);
STD_T = zeros(numel(geneName),1);
STD_C = zeros(numel(geneName),1); 
for i = 1 : numel(geneName)
    STD_T(i) = std(profile(i,strcmp(sampleName,sampleClass(1))));
    STD_C(i) = std(profile(i,strcmp(sampleName,sampleClass(2))));
end

if edgeR
    [FC,pVal] = f.import_edgeR(strcat(dataDir,'edgeR\benchmark_edgeR_mtor.csv'));
else
    % t-test
    FC = meanExpr_T - meanExpr_C;
    FCsort = sort(abs(FC));

    pVal = zeros(numel(FC),1);
    for i = 1 : numel(FC)
        [~,pVal(i)] = ttest2(profile(i,strcmp(sampleName,sampleClass(1))),...
            profile(i,strcmp(sampleName,sampleClass(2))));
    end
end

dispGene = zeros(numel(trackGene),1);
for i = 1 : numel(trackGene)
    tmp = find(strcmp(trackGene(i),geneName),1);
    if ~isempty(tmp)
        dispGene(i) = tmp;
    end
end
%dispGene = sort(dispGene);

sigGene = intersect(find(abs(FC) >= FCcutoff),...
        find(pVal <= pValcutoff));
    
if runTopo
    [outlierDistScore,outlierScore,CER] = generateContours2(profile(:,strcmp(sampleName,sampleClass(1))),...
        profile(:,strcmp(sampleName,sampleClass(2))),gridDensity,numContours,1,targetContainment,removeHighLowExpr);
end
if runMAGE
    [outlierScore,FDR] = MAGE(profile(:,strcmp(sampleName,sampleClass(1))),...
        profile(:,strcmp(sampleName,sampleClass(2))),gridDensity,numContours,1,targetContainment,removeHighLowExpr); 
    
    if ~runTopo
        outlierDistScore = zeros(numel(geneName),1);
    end
end

outlierScoreSort = sort(outlierScore);

sigOutlierGene = find(outlierScore >= outlierScoreCutoff);
sig1D2DGene = intersect(sigGene,sigOutlierGene);

sigGene = sigGene(~ismember(sigGene,sig1D2DGene));
sigOutlierGene = sigOutlierGene(~ismember(sigOutlierGene,sig1D2DGene));
normalGene = ~ismember([1:numel(geneName)],...
    union(sigGene,union(sigOutlierGene,sig1D2DGene)))';


% Set outlierScore to zero for inner genes
%outlierScore(outlierDistScore < outlierDistScoreCutoff) = 0;

figure;
subplot(2,2,1)
scatter(outlierDistScore(normalGene),outlierScore(normalGene),7,[0.7,0.7,0.7],'filled')
hold on
scatter(outlierDistScore(sigGene),outlierScore(sigGene),7,'s','filled')
scatter(outlierDistScore(sigOutlierGene),outlierScore(sigOutlierGene),7,'d','filled')
scatter(outlierDistScore(sig1D2DGene),outlierScore(sig1D2DGene),7,'*')
hold off
yline(outlierScoreCutoff,'--r')
title(strcat('Pearson Correlation: ',num2str(round(corr(outlierDistScore,outlierScore),2))))
xlabel('Dist. from Contour (TPM)')
ylabel('Outlier Score')
subplot(2,2,2)
scatter((STD_T+STD_C)/2,outlierScore,7,[0.7,0.7,0.7],'filled')
hold on
scatter((STD_T(sigGene)+STD_C(sigGene))/2,outlierScore(sigGene),7,'s','filled')
scatter((STD_T(sigOutlierGene)+STD_C(sigOutlierGene))/2,outlierScore(sigOutlierGene),7,'d','filled')
scatter((STD_T(sig1D2DGene)+STD_C(sig1D2DGene))/2,outlierScore(sig1D2DGene),7,'*')
hold off
yline(outlierScoreCutoff,'--r')
title(strcat('Pearson Correlation: ',num2str(round(corr((STD_T+STD_C)/2,outlierScore),2))))
xlabel('Mean STD')
ylabel('Outlier Score')
subplot(2,2,3)
scatter((meanExpr_T+meanExpr_C)/2,outlierScore,7,[0.7,0.7,0.7],'filled')
hold on
scatter((meanExpr_T(sigGene)+meanExpr_C(sigGene))/2,outlierScore(sigGene),7,'s','filled')
scatter((meanExpr_T(sigOutlierGene)+meanExpr_C(sigOutlierGene))/2,outlierScore(sigOutlierGene),7,'d','filled')
scatter((meanExpr_T(sig1D2DGene)+meanExpr_C(sig1D2DGene))/2,outlierScore(sig1D2DGene),7,'*')
hold off
yline(outlierScoreCutoff,'--r')
title(strcat('Pearson Correlation: ',num2str(round(corr((meanExpr_T+meanExpr_C)/2,outlierScore),2))))
xlabel('Mean Expr (TPM)')
ylabel('Outlier Score')

figure;
subplot(2,2,1)
scatter(FC(normalGene),-log(pVal(normalGene)),7,[0.7,0.7,0.7],'filled')
hold on
scatter(FC(sigGene),-log(pVal(sigGene)),7,'s','filled')
scatter(FC(sigOutlierGene),-log(pVal(sigOutlierGene)),7,'d','filled')
scatter(FC(sig1D2DGene),-log(pVal(sig1D2DGene)),7,'*')
yline(-log(pValcutoff),'--r')
if ~ismember(0,dispGene)
    text(FC(dispGene)+0,-log(pVal(dispGene))+0,geneName(dispGene),'FontSize',7)
end
xline(FCcutoff,'--r')
xline(-FCcutoff,'--r')
hold off
title('treatment/control Volcano')
xlabel('log_2(FC)')
ylabel('-log_1_0(p-value)')
legend({'normal','DEG','AEG','DE & AE'})


%figure;
% subplot(2,2,2)
% scatter(find(ismember(sort(FC),FC(normalGene))),sort(FC(normalGene)),9,[0.7,0.7,0.7],'filled')
% hold on
% scatter(find(ismember(sort(FC),FC(sigGene))),sort(FC(sigGene)),15,'s','filled')
% scatter(find(ismember(sort(FC),FC(sigOutlierGene))),sort(FC(sigOutlierGene)),15,'d','filled')
% scatter(find(ismember(sort(FC),FC(sig1D2DGene))),sort(FC(sig1D2DGene)),15,'*')
% yline(FCcutoff,'--r')
% yline(-FCcutoff,'--r')
% hold off
% title('Average FC')
% xticks([])
% xlabel('rank ordered genes')
% ylabel('Log_2(FC)')


meanExpr = mean(profile,2);
%figure;
subplot(2,2,3)
scatter(meanExpr(normalGene),FC(normalGene),9,[0.7,0.7,0.7],'filled')
hold on
scatter(meanExpr(sigGene),FC(sigGene),15,'s','filled')
scatter(meanExpr(sigOutlierGene),FC(sigOutlierGene),15,'d','filled')
scatter(meanExpr(sig1D2DGene),FC(sig1D2DGene),15,'*')
yline(FCcutoff,'--r')
yline(-FCcutoff,'--r')
if ~ismember(0,dispGene)
    text(meanExpr(dispGene),FC(dispGene),geneName(dispGene),'FontSize',7)
end
hold off
title('MA')
xlabel('Mean TPM')
ylabel('log_2(FC)')

X = zeros(numel(outlierScoreSort),1);
for i = 1 : numel(outlierScoreSort)
    if ~isempty(find(ismember(outlierScoreSort(i),outlierScore(normalGene)),1))
        X(i) = i;
    end
end
X = X(X > 0);

%figure;
% subplot(2,2,4)
% scatter(find(X),outlierScoreSort(X),20,[0.7,0.7,0.7],'filled')
% hold on
% scatter(find(ismember(outlierScoreSort,outlierScore(sigGene))),outlierScoreSort(find(ismember(outlierScoreSort,outlierScore(sigGene)))),25,'s','filled')
% scatter(find(ismember(outlierScoreSort,outlierScore(sigOutlierGene))),outlierScoreSort(find(ismember(outlierScoreSort,outlierScore(sigOutlierGene)))),25,'d','filled')
% scatter(find(ismember(outlierScoreSort,outlierScore(sig1D2DGene))),outlierScoreSort(find(ismember(outlierScoreSort,outlierScore(sig1D2DGene)))),25,'*')
% yline(outlierScoreCutoff,'--r')
% %text(find(ismember(outlierScoreSort,outlierScore(dispGene))),outlierScoreSort(find(ismember(outlierScoreSort,outlierScore(dispGene)))),geneName(dispGene),'FontSize',7)
% hold off
% title('Outlier Probability')
% xticks([])
% xlabel('rank ordered genes')
% ylabel('Outlier Probability Score')

%figure;
subplot(2,2,4)
scatter(abs(FC(normalGene)),outlierScore(normalGene),9,[0.7,0.7,0.7],'filled')
hold on
scatter(abs(FC(sigGene)),outlierScore(sigGene),15,'s','filled')
scatter(abs(FC(sigOutlierGene)),outlierScore(sigOutlierGene),15,'d','filled')
scatter(abs(FC(sig1D2DGene)),outlierScore(sig1D2DGene),15,'*')
yline(outlierScoreCutoff,'--k','LineWidth',2)
xline(FCcutoff,'--k','LineWidth',2)
hold off
title(strcat('Spearman Correlation: ',num2str(corr(abs(FC),outlierScore,'Type','Spearman'))))
xlabel('|log_2(FC)|')
ylabel('OS')

figure;
subplot(3,1,1)
histogram(abs(FC))
xline(FCcutoff,'--r','LineWidth',2)
title('Fold-change Distiribution')
xlabel('|log_2(FC)|')
ylabel('Genes')
subplot(3,1,2)
histogram(-log(pVal))
xline(-log(pValcutoff),'--r','LineWidth',2)
title('P-value Distiribution')
xlabel('-log_1_0(p-value)')
ylabel('Genes')
subplot(3,1,3)
histogram(outlierScore)
xline(outlierScoreCutoff,'--r','LineWidth',2)
title('Outlier Score Distiribution')
xlabel('Outlier Score')
ylabel('Genes')


% summary stats
bxgrp = ones(numel(geneName),1);
for i = 1 : numel(geneName)
    if ismember(i,sig1D2DGene)
        bxgrp(i) = 4;
    elseif ismember(i,sigOutlierGene)
        bxgrp(i) = 3;
    elseif ismember(i,sigGene)
        bxgrp(i) = 2;
    end
end
figure;
subplot(1,4,1)
boxplot(outlierScore,bxgrp,'Notch','on','Labels',{'normal','1D','2D','1D/2D'},'PlotStyle','compact')
title('Outlier Score')
ylabel('Outlier Probability')
subplot(1,4,2)
boxplot(abs(FC),bxgrp,'Notch','on','Labels',{'normal','1D','2D','1D/2D'},'PlotStyle','compact')
title('Fold-change')
ylabel('|log_2(FC)|')
subplot(1,4,3)
boxplot(-log(pVal),bxgrp,'Notch','on','Labels',{'normal','1D','2D','1D/2D'},'PlotStyle','compact')
title('P-value')
ylabel('-log_1_0(p-value)')
subplot(1,4,4)
boxplot(log2(meanExpr),bxgrp,'Notch','on','Labels',{'normal','1D','2D','1D/2D'},'PlotStyle','compact')
title('Expression')
ylabel('log_2(TPM)')

% export gene signatures
if exportForPE
    if useTopNgenes > 0
        tmpInd = find(pVal < 0.05);
        [~,tmpSortInd] = sort(abs(FC(tmpInd)));
        sigGenePE = tmpInd(tmpSortInd(end-useTopNgenes+1:end));
        disp("max p-value: "+string(max(abs(pVal(tmpInd(tmpSortInd(end-useTopNgenes+1:end)))))))
        disp("min FC: "+string(min(abs(FC(tmpInd(tmpSortInd(end-useTopNgenes+1:end)))))))
        
        [~,tmpSortInd] = sort(outlierScore);
        sigOutlierGenePE = tmpSortInd(end-useTopNgenes+1:end);
        disp("min OS: "+string(min(outlierScore(tmpSortInd(end-useTopNgenes+1:end)))))
        
        
        sig1D2DGenePE = intersect(sigGenePE,sigOutlierGenePE);
    else
        sigGenePE = sigGene;
        sigOutlierGenePE = sigOutlierGene;
        sig1D2DGenePE = sig1D2DGene;
    end
    
    fid = fopen('geneSig1D','w');
    for i = 1 : numel(sigGenePE)-1
        txt = strcat(geneName(sigGenePE(i)),'\n');
        fprintf(fid,txt);
    end
    txt = [geneName(sigGenePE(i+1))];
    fprintf(fid,txt);
    fclose(fid);
    fid = fopen('geneSig2D','w');
    for i = 1 : numel(sigOutlierGenePE)-1
        txt = strcat(geneName(sigOutlierGenePE(i)),'\n');
        fprintf(fid,txt);
    end
    txt = [geneName(sigOutlierGenePE(i+1))];
    fprintf(fid,txt);
    fclose(fid);
    fid = fopen('geneSig1D2D','w');
    for i = 1 : numel(sig1D2DGenePE)-1
        txt = strcat(geneName(sig1D2DGenePE(i)),'\n');
        fprintf(fid,txt);
    end
    txt = [geneName(sig1D2DGenePE(i+1))];
    fprintf(fid,txt);
    fclose(fid);
    fid = fopen('geneSigBG','w');
    for i = 1 : numel(geneName)-1
        txt = strcat(geneName(i),'\n');
        fprintf(fid,txt);
    end
    txt = [geneName(i+1)];
    fprintf(fid,txt);
    fclose(fid);
    
    clear sigGenePE sigOutlierGenePE sig1D2DGenePE
end

if HC
    sigGeneAll = union(sigGene,union(sigOutlierGene,sig1D2DGene));
    
    clustGeneColors.Labels = cellstr([geneName(sigGene);geneName(sigOutlierGene);geneName(sig1D2DGene)]);
    clustGeneColors.Colors = cellstr([repmat('m',numel(sigGene),1);repmat('b',numel(sigOutlierGene),1);repmat('c',numel(sig1D2DGene),1)]);
    clust = clustergram(profile(sigGeneAll,:),'Symmetric','false',...
        'Colormap','redbluecmap',...
        'OptimalLeafOrder','true',...
        'Displayrange',15,...
        'ColumnLabels',sampleName,...
        'RowLabels',geneName(sigGeneAll),...
        'LabelsWithMarkers','true');
    clust = clustergram([meanExpr_T(sigGeneAll),meanExpr_C(sigGeneAll)],'Symmetric','false',...
        'Colormap','redbluecmap',...
        'OptimalLeafOrder','true',...
        'Displayrange',15,...
        'RowLabels',geneName(sigGeneAll),...
        'ColumnLabels',sampleClass,...
        'LabelsWithMarkers','true');
    
    figure;
    imagesc(profile(sigGeneAll,:))
end

if permFDR
    % permuate samples
    sampleNameOriginal = sampleName;
    for i = 1 : 10
        sampleName = sampleName(randperm(numel(sampleName))');
    end
    
    % Keep only equally expressed (EE) genes
    % by removing genes with highest 5% OS
    FCbelowAlpha = FC(pVal <= 0.05);
    [~,EEind] = sort(abs(FCbelowAlpha));
    EEcutoff = abs(FCbelowAlpha(EEind(round(0.95*numel(EEind)))));
    EEind = find(abs(FC) < EEcutoff | pVal > 0.05);
    
%     % Keep only equally expressed (EE) genes
%     if modAvoidFDROverestimate
%         convEE = find(~ismember(1:size(profile,1),union(sigGene,sig1D2DGene)));
%         topoEE = find(~ismember(1:size(profile,1),union(sigOutlierGene,sig1D2DGene)));
%     else
%         convEE = 1:size(profile,1);
%         topoEE = 1:size(profile,1);
%     end
%     interEE = intersect(convEE,topoEE);
%     convEEint = zeros(numel(interEE),1);
%     topoEEint = zeros(numel(interEE),1);
%     for i = 1 : numel(interEE)
%         convEEint(i) = find(convEE == interEE(i));
%         topoEEint(i) = find(topoEE == interEE(i));
%     end
    
    % permutate expr values within samples
%     for i = 1 : size(profile,2)
%         randInd = randperm(size(profile,1))';
%         profile(:,i) = profile(randInd,i);
%     end
    
    %calculate FC and p using permutated profile
    meanExpr_T_P = mean(profile(EEind,strcmp(sampleName,sampleClass(1))),2);
    meanExpr_C_P = mean(profile(EEind,strcmp(sampleName,sampleClass(2))),2);
    STD_T_P = zeros(numel(EEind),1);
    STD_C_P = zeros(numel(EEind),1);
    for i = 1 : numel(EEind)
        STD_T_P(i) = std(profile(EEind(i),strcmp(sampleName,sampleClass(1))));
        STD_C_P(i) = std(profile(EEind(i),strcmp(sampleName,sampleClass(2))));
    end
    FC_P = meanExpr_T_P - meanExpr_C_P;
    FCsort_P = sort(abs(FC_P));

    pVal_P = zeros(numel(FC_P),1);
    for i = 1 : numel(FC_P)
        [~,pVal_P(i)] = ttest2(profile(EEind(i),strcmp(sampleName,sampleClass(1))),...
            profile(EEind(i),strcmp(sampleName,sampleClass(2))));
    end
    
    %calculate outlier scores using permutated profile
%     [outlierDistScore_P,outlierScore_P] = generateContours2(...
%         profile(topoEE,strcmp(sampleName,sampleClass(1))),...
%         profile(topoEE,strcmp(sampleName,sampleClass(2))),gridDensity,numContours,1);
%     
    %calculate FDR
%     numD = numel(find(abs(FC) >= FCcutoff & pVal <= pValcutoff));
%     numD_P = numel(find(FC_P >= FCcutoff & pVal_P <= pValcutoff));
%     disp(['Permutation test FDR --->   Conventional FDR: ' num2str(100*numD_P/numD) '%'])
%     numD = numel(find(outlierDistScore >= outlierDistScoreCutoff &...
%         outlierScore >= outlierScoreCutoff));
%     numD_P = numel(find(outlierDistScore_P >= outlierDistScoreCutoff &...
%         outlierScore_P >= outlierScoreCutoff));
%     disp(['                            Topological FDR: ' num2str(100*numD_P/numD) '%'])
%     
    figure;
%     subplot(2,2,1)
%     scatter(outlierDistScore_P(topoEEint),outlierScore_P(topoEEint),7,'filled')
%     yline(outlierScoreCutoff,'--r')
%     title(strcat('Pearson Correlation: ',num2str(round(corr(outlierDistScore_P,outlierScore_P),2))))
%     xlabel('Dist. from Contour (TPM)')
%     ylabel('Outlier Score')
%     subplot(2,2,2)
%     scatter((STD_T_P(convEEint)+STD_C_P(convEEint))/2,outlierScore_P(topoEEint),7,'filled')
%     yline(outlierScoreCutoff,'--r')
%     title(strcat('Pearson Correlation: ',num2str(round(corr((STD_T_P(convEEint)+STD_C_P(convEEint))/2,outlierScore_P(topoEEint)),2))))
%     xlabel('Mean STD')
%     ylabel('Outlier Score')
%     subplot(2,2,3)
%     scatter((meanExpr_T_P(convEEint)+meanExpr_C_P(convEEint))/2,outlierScore_P(topoEEint),7,'filled')
%     yline(outlierScoreCutoff,'--r')
%     title(strcat('Pearson Correlation: ',num2str(round(corr((meanExpr_T_P(convEEint)+meanExpr_C_P(convEEint))/2,outlierScore_P(topoEEint)),2))))
%     xlabel('Mean Expr (TPM)')
%     ylabel('Outlier Score')
%     subplot(2,2,4)
    tempCutoff = 0.00001:0.00001:max(abs(FC));
    FDRbyFC = zeros(numel(tempCutoff),2);
    for i = 1 : numel(tempCutoff)
        numD = numel(find(abs(FC(pVal < 0.05)) >= tempCutoff(i)));
        numD_P = numel(find(abs(FC_P(pVal_P < 0.05)) >= tempCutoff(i)));
        FDRbyFC(i,1) = tempCutoff(i);
        FDRbyFC(i,2) = numD_P/(numD_P + numD);
        [i numel(tempCutoff)]
    end
    FDRbyFC(isnan(FDRbyFC(:,2)),2) = 1;
    scatter(tempCutoff',FDRbyFC(:,2),7,'filled')
    title('Sample Permutation FDR')
    xlabel('FC Cutoff')
    ylabel('FDR')
%     sgtitle('Sample Permutation FDR Test')
    
    % assign FDR to genes
        DE_FDR = zeros(size(FC,1),1);
        for i = 1 : numel(DE_FDR)
            [~,tmp] = min(abs(abs(FC(i)) - FDRbyFC(:,1)));
            DE_FDR(i) = FDRbyFC(tmp,2);
        end
    
    
    sampleName = sampleNameOriginal;
end

if saveToTable
    gene_symbol = geneName;
    OS = outlierScore;
    p_value = pVal;
    
    % sort by FDR
    [~,ind] = sort(FDR);
    gene_symbol = gene_symbol(ind);
    OS = OS(ind);
    AE_FDR = FDR(ind);
    FC = FC(ind);
%     p_value = p_value(ind);
    DE_FDR = DE_FDR(ind);
    
    T = table(gene_symbol,OS,AE_FDR,FC,DE_FDR);
    writetable(T,tableName)
end
