
classdef DEG_contour_functions
    
    methods (Static)

%% load gt3
function [profile,geneName,sampleName,grp,grpclr,lgnd,title] = load_gt3()
    % load contour functions
    f = DEG_contour_functions;
    load(strcat(pwd,'\workspaces\breast_cancer_abemaciclib_NCBI.mat'));

    geneName = geneName_GDS4059fullabemaciclib;
    sampleName = sampleName_GDS4059fullabemaciclib;
    profile = profile_GDS4059fullabemaciclib;
    clear geneName_GDS4059fullabemaciclib sampleName_GDS4059fullabemaciclib profile_GDS4059fullabemaciclib

    % remove last gene which is NaN
    profile = profile(1:end-1,:);
    geneName = geneName(1:end-1);

    % combine duplicate genes (microarray probes)
    geneName_all = geneName;
    profile_all = profile;
    geneName = unique(geneName);
    profile = zeros(numel(geneName),size(profile_all,2));
    for i = 1 : numel(geneName)
        ind = find(strcmp(geneName(i),geneName_all));
        profile(i,:) = mean(profile_all(ind,:));
    end
    clear profile_all geneName_all ind

    % remove non genecode genes
    geneName = geneName(416:end);
    profile = profile(416:end,:);

    % create legend, group, and colors
    lgnd = unique(sampleName);

    grp = zeros(numel(sampleName),1);
    for i = 1 : numel(sampleName)
        grp(i) = find(strcmp(sampleName(i),lgnd));
    end
    clrset = hot(numel(lgnd));
    grpclr = clrset;

    title = "breast cancer (abemaciclib)";

    %f.filterAnalysis(profile,[1,3,5,8],[1,2,3,6]);
    [~, ind] = f.fltr(profile,1,6);
    profile = profile(ind,:);
    geneName = geneName(ind);

    % log shift
    %profile = f.logNorm(profile,1,0);
end

%% load mTOR
function [profile,geneName,sampleName,grp,grpclr,lgnd,title] = load_mTOR()
    % load contour functions
    f = DEG_contour_functions;
    
    load(strcat(pwd,'\workspaces\breast_cancer_mouse_mTOR_delta_NCBI.mat'));

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

    title = "breast cancer (mouse mTOR KO)";

    %[final_profile, filteredInd] = fltr(profile,level,sampleMin)
    %f.filterAnalysis(profile,[0.1,1,3,5,8],[1,2,3,6]);
    [~, ind] = f.fltr(profile,0.1,3);
    profile = profile(ind,:);
    geneName = geneName(ind);
end
        
%% benchmark test: Agreement (11/9/2023)
function benchTestAgreement(test1,test2,Sweepstep)
    % Overlap w/ test parameter sweep
    
    %             intersect (test_1 < C) & (test_2 < C)
    % overlap =   -------------------------------------
    %                       min (# of genes)
    
    cutoff = [0 : Sweepstep : 1];
    overlap = zeros(numel(cutoff),numel(cutoff));
    for i = 1 : numel(cutoff)
        test1_positives = find(test1 < cutoff(i));
        for j = 1 : numel(cutoff)
            test2_positives = find(test2 < cutoff(j));
            overlap(i,j) = numel(intersect(test1_positives,test2_positives))/...
                min(numel(test1_positives),numel(test2_positives));
        end
    end
    
    figure;
    s = surf(overlap*100);
    axis([1 size(overlap,1) 1 size(overlap,2)])
    xticklabels(round(str2double(xticklabels)./size(overlap,1),2))
    yticklabels(round(str2double(yticklabels)./size(overlap,2),2))
    if numel(cutoff) > 300
        s.EdgeColor = 'none';
    end
end

%% import edgeR results
function [FC,pVal] = import_edgeR(fileName)
    % import data from R (.csv) (11/6/2023)
    % After profile has been loaded workspace should contain:
    % profile ---> contains un-normalized RNA-seq read counts
    % geneName ---> list of gene names
    % sampleName ---> list of sample names

    benchmark_edgeR = readmatrix(fileName,'OutputType','string');

    geneName_edgeR = benchmark_edgeR(:,1);
    FC = str2double(benchmark_edgeR(:,2));
    CPM = str2double(benchmark_edgeR(:,3));
    Lratio = str2double(benchmark_edgeR(:,4));
    pVal = str2double(benchmark_edgeR(:,5));
    F = str2double(benchmark_edgeR(:,6));
    FpVal = str2double(benchmark_edgeR(:,7));
    clear benchmark_edgeR

    figure;
    subplot(2,2,1)
    scatter(FC,-log10(pVal))
    xlabel('log_2(FC)')
    ylabel('-log_1_0(p-value)')
    subplot(2,2,2)
    scatter(CPM,FC)
    xlabel('CPM')
    ylabel('log_2(FC)')
    subplot(2,2,3)
    scatter(Lratio,F)
    xlabel('likelihood - ratio')
    ylabel('F - test')
    subplot(2,2,4)
    scatter(-log10(pVal),-log10(FpVal))
    xline(-log10(0.05),'--r')
    yline(-log10(0.05),'--r')
    xlabel('-log_1_0(p-value)')
    ylabel('-log_1_0(quasi-likelihood p-value)')
    sgtitle('benchmark - edgeR - output')
    
end

%% filter analysis
function filterAnalysis(profile,level,sampleMin)
% filter: remove genes w/ less than 'filter' # samples with > 1 read
%
% log: log2(X + 1) shift if log > 0
%
% norm: 0) none
%       1) normalize samples
%       2) normalize genes
%       3) normalize samples & genes
%
%
%
    keepGene = false(size(profile,1),1);
    numGeneKept = zeros(numel(level),numel(sampleMin));
    for l = 1 : numel(level)
        for s = 1 : numel(sampleMin)
            
            for g = 1 : size(profile,1)
                keepGene(g) =...
                    numel(find(profile(g,:) > level(l))) >= sampleMin(s);
                [l s g] %#ok<NOPRT>
            end
            numGeneKept(l,s) = numel(find(keepGene));
            
        end
    end
    
    figure;
    for i = 1 : size(numGeneKept,1)
        plot(sampleMin,numGeneKept(i,:),'LineWidth',2)
        hold on
    end
    hold off
    title('filter parameters')
    xlabel('min number of samples > level')
    ylabel('number of genes passing filter')
    legend(strcat('lvl: ',num2str(level')))
end

%% filter low expr
function [final_profile, filteredInd] = fltr(profile,level,sampleMin) 
    keepGene = false(size(profile,1),1);
    for g = 1 : size(profile,1)
        keepGene(g) =...
            numel(find(profile(g,:) > level)) >= sampleMin;
    end
    final_profile = profile(keepGene,:);
    filteredInd = keepGene;
end
        
%% log shift --> normalize

function final_profile = logNorm(profile,log,norm)
%
% log: log2(X + 1) shift if log > 0
%
% norm: 0) none
%       1) normalize samples
%       2) normalize genes
%       3) normalize samples & genes
%
%
%



    % log2 shift of read counts
    if log > 0
        profile = log2(profile + 1);

        figure;
        subplot(2,1,1)
        histogram(mean(profile,1))
        title('expression per sample')
        xlabel('mean log2(read counts + 1)')
        ylabel('sample count')
        subplot(2,1,2)
        histogram(mean(profile,2))
        %xline(cutoff,'-r',{'low read limit'})
        title('expression per gene')
        xlabel('mean log2(read counts + 1)')
        ylabel('gene count')
        set(gca,'Yscale','log')
        sgtitle('log shifted read counts')
    end

    % normalize each cell
    if norm == 1 || norm == 3
        profile = normalize(profile,1);

        figure;
        subplot(2,1,1)
        histogram(mean(profile,1))
        title('expression per sample')
        xlabel('mean log2(read counts + 1)')
        ylabel('sample count')
        subplot(2,1,2)
        histogram(mean(profile,2))
        title('expression per gene')
        xlabel('mean log2(read counts + 1)')
        ylabel('gene count')
        set(gca,'Yscale','log')
        sgtitle('normalized across samples')
    end

    % normalize each gene
    if norm == 2 || norm == 3
        profile = normalize(profile,2);

        figure;
        subplot(2,1,1)
        histogram(mean(profile,1))
        title('expression per sample')
        xlabel('mean log2(read counts + 1)')
        ylabel('sample count')
        subplot(2,1,2)
        histogram(mean(profile,2))
        title('expression per gene')
        xlabel('mean log2(read counts + 1)')
        ylabel('gene count')
        set(gca,'Yscale','log')
        sgtitle('normalized across genes')
    end
    
    final_profile = profile;
end

%% Function: sim ellipse

function [] = simEllipse(data1, data2, ind, lowercutoff, uppercutoff,...
    gridDensity, saveGenes, geneids)

    corr_matrix = GeneSimilarityMatrix(data1(ind,:)',data2(ind,:)');

    tissueind = selectPairs(corr_matrix, lowercutoff, uppercutoff);

    timesOut = outlierDist(data1,data2,tissueind,gridDensity,1);

    if saveGenes
        saveOutliers(timesOut, geneids);
    end
    
    % output info
    size(tissueind,1)
end

%% Function: sim contour

function [] = simContour(data1, data2, ind, lowercutoff, uppercutoff,...
    gridDensity, saveGenes, geneids)

    corr_matrix = GeneSimilarityMatrix(data1(ind,:)',data2(ind,:)');

    tissueind = selectPairs(corr_matrix, lowercutoff, uppercutoff);

    timesOut = outlierDist(data1,data2,tissueind,gridDensity,0);

    if saveGenes
        saveOutliers(timesOut, geneids);
    end
    
    % output info
    size(tissueind,1)
end

%% Function: sim residual

function [] = simResidual(data1, data2, ind, lowercutoff, uppercutoff,...
    saveGenes, geneids)

    corr_matrix = GeneSimilarityMatrix(data1(ind,:)',data2(ind,:)');

    tissueind = selectPairs(corr_matrix, lowercutoff, uppercutoff);
    
    geneResiduals = zeros(numel(ind),1);
    
    for i = 1 : size(tissueind,1)
        geneResiduals = geneResiduals + abs(residualAnalysis(...
            data1(ind,tissueind(i,1)), data2(ind,tissueind(i,1)),0));
    end
    
    figure;
    histogram(geneResiduals,'FaceColor','g');
    title('Distribution of residual sums');
    xlabel('Residual sum');
    ylabel('# of genes');


    if saveGenes
        saveOutliers(geneResiduals, geneids);
    end
    
    % output info
    size(tissueind,1)
end

%% Function: bin ellipse/contour and residual (8/15/2022)

function [] = binER(data1, data2, ind, lowercutoff, uppercutoff,...
    gridDensity, geneids, makeEllipse)

    corr_matrix = GeneSimilarityMatrix(data1(ind,:)',data2(ind,:)');

    tissueind = selectPairs(corr_matrix, lowercutoff, uppercutoff);
    
    geneResiduals = zeros(numel(ind),1);
    
    for i = 1 : size(tissueind,1)
        geneResiduals = geneResiduals + abs(residualAnalysis(...
            data1(ind,tissueind(i,1)), data2(ind,tissueind(i,1)),0));
    end
    
    timesOut = outlierDist(data1,data2,tissueind,gridDensity,makeEllipse);
    
    [~,sortedOutlierGeneIndE] = sort(timesOut,1,'descend');
    outlierGeneIndE = sortedOutlierGeneIndE(1:round(numel(timesOut)*0.05));
    [~,sortedOutlierGeneIndR] = sort(geneResiduals,1,'descend');
    outlierGeneIndR = sortedOutlierGeneIndR(1:round(numel(geneResiduals)*0.05));
    
    % bin
    outlierBin = ones(numel(geneids),1);
    outlierER = intersect(outlierGeneIndE,outlierGeneIndR);
    outlierBin(outlierER) = 4;
    outlierBin(outlierGeneIndE(~ismember(outlierGeneIndE,outlierER))) = 2;
    outlierBin(outlierGeneIndR(~ismember(outlierGeneIndR,outlierER))) = 3;
    
    % save gene Bins
    for bin = 1 : 4
        % file name
        switch bin
            case 1
                fid = fopen('outlierGenesNone.txt','w');
            case 2
                fid = fopen('outlierGenesE.txt','w');
            case 3
                fid = fopen('outlierGenesR.txt','w');
            case 4
                fid = fopen('outlierGenesER.txt','w');
        end
        
        binInd = find(outlierBin == bin);
        
        for i = 1 : numel(binInd)-1
            text = strcat(geneids(binInd(i)),'\n');
            fprintf(fid,text);
        end
        text = [geneids(binInd(i+1))];
        fprintf(fid,text);
        fclose(fid);
    end
end

%% Function: select highly correlated/anticorrelated tissue/line pairs (6/23/2021)

function tissueind = selectPairs(corr_matrix,lowercutoff, uppercutoff)


tissueind = zeros(1,2);
for i = 1 : size(corr_matrix,1)
    for j = 1 : size(corr_matrix,2)
        if abs(corr_matrix(i,j)) > lowercutoff &&...
                abs(corr_matrix(i,j)) < uppercutoff
            tissueind = [tissueind ; i,j];
        end
    end
end
tissueind = tissueind(2:end,:);

end

%% Function: outlier distribution
function timesOut = outlierDist(data1, data2,repGrp1, repGrp2, gridDensity, method, disp)
    
    if disp
        figure;
    end
    
    if method == 1
        [timesOut,~] =...
            generateEllipses(data1,data2,gridDensity);
    elseif method == 2
        [timesOut,~] =...
            generateContours(data1,data2,repGrp1,gridDensity,disp);
    elseif method == 3
        [timesOut,~] =...
            generateContours2(data1,data2,gridDensity,disp);
    end
    
    if disp
        figure;
        n = max(timesOut);
        histogram(timesOut,n+1,'FaceColor','g');
        title('Distribution of outliers');
        xlabel('Occurances outside of normal region');
        ylabel('# of genes');

        if n < 100
            xticks(n/(n+1)/2:n/(n+1):n-n/(n+1)/2);
            xticklabels(0:n);
        end
    end

end

%% Function: save outlier gene names (8/3/2022)
function [] = saveOutliers(timesOut, geneid)

    [~,sortedOutlierGeneInd] = sort(timesOut,1,'descend');
    outlierGeneInd = sortedOutlierGeneInd(1:round(numel(timesOut)*0.05));

    fid = fopen('outlierGenes.txt','w');
    for i = 1 : numel(outlierGeneInd)-1
        text = strcat(geneid(outlierGeneInd(i)),'\n');
        fprintf(fid,text);
    end
    text = [geneid(outlierGeneInd(i+1))];
    fprintf(fid,text);
    fclose(fid);

end

%% Function: display ellipse method (8/9/2022)
% Create contour plot as easier method for creating ellipse (6/28/2021)

function [] = displayEllipseStep(gtex,ccle,i,j, gridDensity)

    % display points
    figure;
    subplot(2,3,1);
    scatter(ccle(:,j),gtex(:,i),25,'r.');
    axis equal
    title('Selected tissue/cell line pair');
    %xlabel(newccle_tissue{1,j});
    %ylabel(newgtex_tissue{1,i});

    % create a grid
    gridheight = range(gtex(:,i))/gridDensity;
    gridwidth = range(ccle(:,j))/gridDensity;
    densityMat = zeros(gridDensity,gridDensity);
    startx = min(ccle(:,j)) + gridwidth + 0.00001; % correction
    starty = min(gtex(:,i)) + gridheight + 0.00001;

    % find density within each grid region
    for k = 1 : size(gtex,1)
        % start at 1,1 in grid
        gridx = startx;
        gridy = starty;
        addAtx = 1;
        addAty = 1;

        % find where point is in grid
        while 1
            if ccle(k,j) < gridx
                while 1
                    if gtex(k,i) < gridy
                        break
                    end
                    gridy = gridy + gridheight;
                    addAty = addAty + 1;
                end
                break
            end
            gridx = gridx + gridwidth;
            addAtx = addAtx + 1;
        end

        % add point to density matrix
        densityMat(addAty,addAtx) = densityMat(addAty,addAtx) + 1;
    end

    % display density matrix
    subplot(2,3,2);
    imagesc(flip(densityMat,1));
    title('Density Matrix');
    axis square
    xticks([]);
    yticks([]);
    ylabel('arb. units');

    % make contour plot of density matrix
    subplot(2,3,3);
    [densityContPoints,densityCont] = contour(densityMat,20);
    axis equal
    title('Contours');
    xticks([]);
    yticks([]);
    ylabel('arb. units');

    % find 95% confidence contour
    lvls = densityCont.LevelList;
    confidenceContInd = find(densityContPoints(1,:) == lvls(2)); % indices at contour level 2 ~ approx. 95%
    numVertices = densityContPoints(2,confidenceContInd);
    confidenceCont = zeros(numVertices,2);
    for k = 1 : numVertices
        confidenceCont(k,1) = densityContPoints(1,confidenceContInd+k);
        confidenceCont(k,2) = densityContPoints(2,confidenceContInd+k);
    end

    subplot(2,3,4);
    plot(confidenceCont(:,1),confidenceCont(:,2));
    hold on
    title('95% containment contour');

    % create an ellipse based on 95% contour

    % find major axis along widest direction
    point1 = [confidenceCont(1,1) confidenceCont(1,2)];
    point2 = point1;
    maxPoints = [point1 point2];
    majorLength = 0;
    for k = 1 : size(confidenceCont,1)
        for k2 = 1 : size(confidenceCont,1)
            dist = pdist([point1;point2]);
            if dist > majorLength
                majorLength = dist;
                maxPoints = [point1;point2];
            end
            point2 = [confidenceCont(k2,1) confidenceCont(k2,2)];
        end
        point1 = [confidenceCont(k,1) confidenceCont(k,2)];
    end
    majorLength = majorLength/2;
    plot(maxPoints(:,1),maxPoints(:,2),'g');


    % ellipse center is midpoint of major axis
    cent = [(maxPoints(2,1)/2 + maxPoints(1,1)/2),(maxPoints(2,2)/2 + maxPoints(1,2)/2)];
    s = [1;(maxPoints(2,2)-maxPoints(1,2))/(maxPoints(2,1)-maxPoints(1,1))];
    yIntercept = maxPoints(1,2) - s(2,1)*maxPoints(1,1); % y intercept: b = y - mx 
    plot(cent(1,1), cent(1,2),'x','MarkerSize',10,'LineWidth',2,'Color','r');

    % determine azmuthal angle from slope
    s = s(2,1);
    azmuth = 180*atan(s)/pi;  % convert slope to azimuthal angle

    % find minor axis length
    minorLength = abs(p_poly_dist(cent(1),cent(2),confidenceCont(:,1),confidenceCont(:,2)));
    minorPoints = [cent(1) + minorLength*cos(pi/2+atan(s)) cent(2) + minorLength*sin(pi/2-atan(s));...
                   cent(1) + minorLength*cos(pi/2-atan(s)) cent(2) - minorLength*sin(pi/2+atan(s))];
    plot(minorPoints(:,1),minorPoints(:,2),'g');

    correctRatio =  1;
    while correctRatio

        % create ellipse
        ecc = sqrt(1-((minorLength^2)/(majorLength^2)));  % eccentricity
        [xEarb,yEarb] = ellipse1(cent(1),cent(2),[majorLength ecc],azmuth,[]);  % generate coordinate points of ellipse



        % convert points back to original scale
        xE = xEarb(:)*gridwidth + startx - gridwidth;
        yE = yEarb(:)*gridheight + starty - gridheight;

        in = inpolygon(ccle(:,j),gtex(:,i),xE,yE); % color points in/out of ellipse
        out = find(in(:) == 0);
        in = find(in);

        % make ellipse larger (to encompase entire contour)
        if numel(in)/(numel(in)+numel(out)) >= 0.95
            correctRatio = 0;
        else
            minorLength = minorLength + minorLength*0.1;
            majorLength = majorLength + majorLength*0.1;
        end
    end
    plot(xEarb,yEarb,'g','LineWidth',3);
    axis equal
    xticks([]);
    yticks([]);
    ylabel('arb. units');
    hold off

    subplot(2,3,5);
    scatter(ccle(out,j),gtex(out,i),25,'r.');
    hold on
    scatter(ccle(in,j),gtex(in,i),25,'c.');
    plot(xE,yE,'g','LineWidth',3);
    hold off
    axis equal
    title('Selected tissue/cell line pair w/ confidence ellipse');
    %xlabel(newccle_tissue{1,j});
    %ylabel(newgtex_tissue{1,i});

end

%% Function: display contour method (8/16/2022)
% Create contour plot as easier method for creating ellipse (6/28/2021)

function [] = displayContourStep(gtex,ccle,i,j, gridDensity)

    % display points
    figure;
    tiledlayout(2,3);
    nexttile;
    scatter(ccle(:,j),gtex(:,i),25,'r.');
    axis equal
    title('Selected tissue/cell line pair');
    %xlabel(newccle_tissue{1,j});
    %ylabel(newgtex_tissue{1,i});

    % create a grid
    gridheight = range(gtex(:,i))/gridDensity;
    gridwidth = range(ccle(:,j))/gridDensity;
    densityMat = zeros(gridDensity,gridDensity);
    startx = min(ccle(:,j)) + gridwidth + 0.00001; % correction
    starty = min(gtex(:,i)) + gridheight + 0.00001;

    % find density within each grid region
    for k = 1 : size(gtex,1)
        % start at 1,1 in grid
        gridx = startx;
        gridy = starty;
        addAtx = 1;
        addAty = 1;

        % find where point is in grid
        while 1
            if ccle(k,j) < gridx
                while 1
                    if gtex(k,i) < gridy
                        break
                    end
                    gridy = gridy + gridheight;
                    addAty = addAty + 1;
                end
                break
            end
            gridx = gridx + gridwidth;
            addAtx = addAtx + 1;
        end

        % add point to density matrix
        densityMat(addAty,addAtx) = densityMat(addAty,addAtx) + 1;
    end

    % display density matrix
    nexttile;
    imagesc(flip(densityMat,1));
    title('Density Matrix');
    axis square
    xticks([]);
    yticks([]);
    ylabel('arb. units');

    % make contour plot of density matrix
    nexttile;
    [densityContPoints,densityCont] = contour(densityMat,40);
    axis equal
    title('Contours');
    xticks([]);
    yticks([]);
    ylabel('arb. units');
    
    lvls = densityCont.LevelList;
    % indices at contour level 1 ~ approx. 95%
    confidenceContInd = find(densityContPoints(1,:) == lvls(1));
    numVertices = densityContPoints(2,confidenceContInd);
    ax4 = nexttile;
    
    for region = 1 : numel(numVertices)
        
        confidenceCont = zeros(numVertices(region),2);
        for k = 1 : numVertices(region)
            confidenceCont(k,1) = densityContPoints(1,confidenceContInd(region)+k);
            confidenceCont(k,2) = densityContPoints(2,confidenceContInd(region)+k);
        end

        plot(ax4,confidenceCont(:,1),confidenceCont(:,2),'k');

        
        


         % convert points back to original scale
         confidenceCont(:,1) = confidenceCont(:,1)*gridwidth + startx - gridwidth;
         confidenceCont(:,2) = confidenceCont(:,2)*gridheight + starty - gridheight;

         in = inpolygon(ccle(:,j),gtex(:,i),...
             confidenceCont(:,1),confidenceCont(:,2)); % color points in/out of contour
         out = find(in(:) == 0);
         in = find(in);
         if region == 1
             allin = in;
             hold on
             title(ax4,'95% containment contour');
         else
             allin = union(allin,in);
         end
    end
    numel(allin)/(numel(in)+numel(out))
    
    ax5 = nexttile;
    scatter(ax5,ccle(out,j),gtex(out,i),25,'r.');
    hold on
    scatter(ccle(allin,j),gtex(allin,i),25,'c.');
    plot(confidenceCont(:,1),confidenceCont(:,2),'g','LineWidth',3);
    hold off
    axis equal
    title('Selected tissue/cell line pair w/ confidence ellipse');
    %xlabel(newccle_tissue{1,j});
    %ylabel(newgtex_tissue{1,i});

end

%% Function: residual analysis (8/10/2022)
function geneResiduals = residualAnalysis(data1,data2,dispSteps)

    linReg = fitlm(data2,data1);
    
    if dispSteps
        figure;
        subplot(2,2,1)
        plot(linReg)
        title('Tissue/cell line pair expression')
        ylabel('Tissue i')
        xlabel('Cell line j')
        subplot(2,2,2)
        plotResiduals(linReg,'caseorder');
        subplot(2,3,4)
        plotResiduals(linReg)
        subplot(2,3,5)
        plotResiduals(linReg,'probability')
        subplot(2,3,6)
        plotResiduals(linReg,'symmetry')
        sgtitle('Linear Fit Model')
    end
    
    geneResiduals = linReg.Residuals.Raw;
end

%% Function: oncogene comparison (8/17/2022)

function [pVal,ratioSlope] = oncComp(timesOut,GeneNames,dataDir, disp)

    % find oncogene ind
    load(strcat(dataDir,'cosmic_oncogene_data'))
    oncNames = cosmicCancerGeneSymbol;

    oncInd = zeros(numel(oncNames),1);
    for i = 1 : numel(oncInd)
        temp = find(strcmp(oncNames(i),GeneNames));
        if ~isempty(temp)
            oncInd(i) = temp;
        end
    end
    oncInd = oncInd(oncInd ~= 0);
    nonOncInd = ~ismember(1:numel(GeneNames),oncInd)';

    n = 10; %max(timesOut)+1;

    if disp
        figure;
        subplot(2,2,1)
        hOnc = histogram(timesOut(oncInd),n,'FaceColor','r');
        subplot(2,2,2)
        hnonOnc = histogram(timesOut(nonOncInd),n,'FaceColor','g');
        hOncVal = hOnc.Values;
        hnonOncVal = hnonOnc.Values;
        bx = zeros(numel(timesOut),2);
        bx(1:numel(timesOut(oncInd)),:) = [timesOut(oncInd),...
            ones(numel(timesOut(oncInd)),1)];
        bx(numel(timesOut(oncInd))+1:end,:) = [timesOut(nonOncInd),...
            2*ones(numel(timesOut(nonOncInd)),1)];
    else
        [hOncVal,~] = histcounts(timesOut(oncInd),n);
        [hnonOncVal,~] = histcounts(timesOut(nonOncInd),n);
    end
        
        oncRatio = hOncVal;
        for i = 1 : numel(oncRatio)
            oncRatio(i) = oncRatio(i)/(hOncVal(i) + hnonOncVal(i));
        end
        
    if disp
        subplot(2,1,1)
        boxplot(bx(:,1),bx(:,2),'Orientation','horizontal','Notch','on');
        title('Outlier distribution');
        yticklabels({'Onc','non-onc'});

        subplot(2,1,2)
        scatter(0:n-1,oncRatio,'filled','k')
        title('Oncogene ratio per bin');
        xlabel('Occurances outside of normal region');
        ylabel('# of onc / bin total');
    end

    % p-value
    [~,pVal] = ttest2(timesOut(oncInd),timesOut(nonOncInd));
    
    n = find(isfinite(oncRatio));
    oncRatio = oncRatio(n);
    % ratioSlope
    ratioSlope = polyfit(n-1,oncRatio,1);
    ratioSlope = ratioSlope(1);

end

%% Function: find matching tissues (9/27/2022)
% find members of list1 and 2 which have 4 character match

function matchInd = matchTissues(list1,list2)
    matchInd = zeros(numel(list1)*numel(list2),2);
    temp = 1;
    for i = 1 : numel(list1)
        for j = 1 : numel(list2)
            
            if strlength(list1(i)) >= 4
                k = 1;
                while k <= strlength(list1(i)) - 3
                    pat = extractBetween(list1(i),k,k+3);
                    if contains(list2(j),pat,'IgnoreCase',true)
                        matchInd(temp,1) = i;
                        matchInd(temp,2) = j;
                        temp = temp + 1;
                        k = strlength(list1(i));
                    end
                    k = k + 1;
                end
            end
        end
    end
    matchInd = matchInd(find(matchInd(:,1)),:);
end

%% Function: find matching genes (4/15/2024)
% find members of list1 and 2 which have complete match

function matchInd = matchGenes(list1,list2)
    list1 = lower(list1);
    list2 = lower(list2);
    
    matchInd = zeros(numel(list1)*numel(list2),2);
    temp = 1;
    for i = 1 : numel(list1)
        for j = 1 : numel(list2)
            if strcmp(list1(i),list2(j))
                matchInd(temp,1) = i;
                matchInd(temp,2) = j;
                temp = temp + 1;
            end
        end
    end
    matchInd = matchInd(find(matchInd(:,1)),:);
end

    end
end