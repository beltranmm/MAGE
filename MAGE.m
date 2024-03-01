% generate confidence contours based on (ind) given


function [OutlierDistScore, OutlierScore, printConfidenceCont] = generateContours2(dataX,dataY,varargin)    
    narginchk(2,10)
    if nargin >= 3
        gridDensity = varargin{1};
    else
        gridDensity = 50;
    end
    if nargin >= 4
        numContours = varargin{2};
    else
        numContours = 20;
    end
    if nargin >= 5
        outputPlots = varargin{3};
    else
        outputPlots = false;
    end
    if nargin >= 6
        targetContainment = varargin{4};
    else
        targetContainment = 0.95;
    end
    if nargin >= 7
        removeHighLowExpr = varargin{5};
    else
        removeHighLowExpr = true;
    end
    if nargin >= 8
        contourLoopsMax = varargin{6};
    else
        contourLoopsMax = 10;
    end
    if nargin >= 9
        numStartingContours = varargin{7};
    else
        numStartingContours = 200;
    end
    if nargin >= 10
        distributionSelection = varargin{8};
    else
        distributionSelection = 'gauss';
    end
    
    % validate robust contours
    validation = 0;
    validateAllContours = 0;
    validateN = 10;
    validationOScutoff = 0.5;
    

    numGene = size(dataY,1);
    nX = size(dataX,2);
    nY = size(dataY,2);
    
    contGridAreaRatioCutoff = 0.5;
    geneProbContainedCutoff = 0.7;
    
    disp(['Calculating ' char(956) ' and ' char(963) ' per gene '])
    
    geneMeanX = zeros(numGene,1);
    geneMeanY = zeros(numGene,1);
    geneSTDX = zeros(numGene,1);
    geneSTDY = zeros(numGene,1);
        
    %calculate gene mean and STD across replicate groups
    for k = 1 : numGene
        geneMeanX(k) = mean(dataX(k,:));
        geneMeanY(k) = mean(dataY(k,:));
        
        geneSTDX(k) = std(dataX(k,:)) + 0.0001;
        geneSTDY(k) = std(dataY(k,:)) + 0.0001;
    end
    
    totalSTD = mean([geneSTDX,geneSTDY]);
    
    disp('Fitting a grid ')
    
    if outputPlots
        % display points
        figure;
        tiledlayout(2,3);
        nexttile;
        scatter(geneMeanX,geneMeanY,25,'r.');
        axis equal
        title('Normalized Expression');
        xlabel('TPM');
        ylabel('TPM');
    end


    % create a grid
    gridheight = (range(geneMeanY) + 2*2*mean(abs(geneSTDY)))/gridDensity;
    gridwidth = (range(geneMeanX) + 2*2*mean(abs(geneSTDX)))/gridDensity;
    densityMat = zeros(gridDensity,gridDensity);
    startx = min(geneMeanX) - 2*mean(abs(geneSTDX)) + 0.00001;
    starty = min(geneMeanY) - 2*mean(abs(geneSTDY)) + 0.00001;
    
    disp('Summing gene probabilities')
    % find density within each grid region
    for k = 1 : numGene
        % start at 1,1 in grid
        gridLocationX = startx;
        for gridIndX = 1 : gridDensity

            gridLocationY = starty;
            for gridIndY = 1 : gridDensity
                if strcmp(distributionSelection, 'poisson')
                    densityMat(gridIndY,gridIndX) =...
                        densityMat(gridIndY,gridIndX) +...
                        (poisspdf(round(gridLocationX),round(geneMeanX(k))) *...
                        poisspdf(round(gridLocationY),round(geneMeanY(k))));
                elseif strcmp(distributionSelection, 'gamma')
                    betaX = geneMeanX(k)/(max(geneSTDX(k),1)^2);
                    betaY = geneMeanY(k)/(max(geneSTDX(k),1)^2);
                    alphaX = (geneMeanX(k)/max(geneSTDX(k),1))^2;
                    alphaY = (geneMeanY(k)/max(geneSTDX(k),1))^2;
                    densityMat(gridIndY,gridIndX) =...
                        densityMat(gridIndY,gridIndX) +...
                        (gampdf(gridLocationX,alphaX,betaX) *...
                        gampdf(gridLocationY,alphaY,betaY));
                elseif contains(distributionSelection, 'gaussTrunc')
                    sC = str2double(extractAfter(distributionSelection,'Trunc'));
                    if gridLocationX <= geneMeanX(k)+sC*geneSTDX(k) &&...
                            gridLocationY <= geneMeanY(k)+sC*geneSTDY(k)
                        densityMat(gridIndY,gridIndX) =...
                            densityMat(gridIndY,gridIndX) +...
                            (normpdf(gridLocationX,geneMeanX(k),geneSTDX(k)) *...
                            normpdf(gridLocationY,geneMeanY(k),geneSTDY(k)));
                    end
                else
                    densityMat(gridIndY,gridIndX) =...
                        densityMat(gridIndY,gridIndX) +...
                        (normpdf(gridLocationX,geneMeanX(k),geneSTDX(k)) *...
                        normpdf(gridLocationY,geneMeanY(k),geneSTDY(k)));
                end
                if isnan(densityMat(gridIndY,gridIndX))
                    [k gridIndX gridIndY]
                    break;
                end

                gridLocationY = gridLocationY + gridheight;
            end
            gridLocationX = gridLocationX + gridwidth;
        end
    end
    %densityMat = log2(densityMat+1);


    % make contour plot of density matrix
    if outputPlots
        % display density matrix
        nexttile;
        imagesc(flip(densityMat,1));
        title('Density Matrix');
        axis square
        xticks([]);
        yticks([]);

        ax3 = nexttile;
    end

    % pad boundary to ensure closed contour lines
    disp('Selecting characteristic contour ')
    densityMat = padarray(densityMat,[1 1], mean(densityMat,"all"), 'both');
    contourLoopsOptimalInfo = zeros(contourLoopsMax,3);
    contourRange = numStartingContours;
    countourCutoff = round(0.25*numContours);
    c = 0; % track number of contour loops
    
    while c < contourLoopsMax
        c = c + 1;
        if outputPlots
            [densityContPoints,densityCont] = ...
                contour(ax3,log(densityMat+1),contourRange);
            hold(ax3,'on')
        else
            [densityContPoints,densityCont] = ...
                contour(log(densityMat+1),contourRange);
        end
        
        % evaluate contour effectivness
        lv = 0;
        lvls = densityCont.LevelList;

        contourIsNotCorrect = 1;
        zeroVertices = 1;
        testingContourSelection = 1;
        geneProbContainedInLv = zeros(numel(lvls),1);
        areaRatioOfLv = zeros(numel(lvls),1);

        while contourIsNotCorrect || zeroVertices || testingContourSelection
            lv = lv + 1;
            confidenceContInd = zeros(size(densityContPoints,2),1);
            temp = 1;
            while temp < size(densityContPoints,2)
                if densityContPoints(1,temp) == lvls(lv)
                    confidenceContInd(temp) = 1;
                end
                temp = temp + densityContPoints(2,temp) + 1;
            end

            confidenceContInd = find(confidenceContInd);
            numVertices = densityContPoints(2,confidenceContInd);

            if numel(numVertices) > 0

                zeroVertices = 0;

                % contour:grid area ratio < contGridAreaRatioCutoff
                monteCarloSampleSize = 1000;
                randPoints = [gridDensity*rand(monteCarloSampleSize,1),...
                    gridDensity*rand(monteCarloSampleSize,1)];
                monteCarloSmallSampleSize = 10;
                randGenePoints = zeros(numGene*monteCarloSmallSampleSize,2);
                temp = 1;
                for k = 1 : numGene
                    randGenePoints(temp:temp+monteCarloSmallSampleSize-1,:) =...
                        [geneMeanX(k) + geneSTDX(k)*randn(monteCarloSmallSampleSize,1),...
                        geneMeanY(k) + geneSTDY(k)*randn(monteCarloSmallSampleSize,1)];
                    temp = temp + monteCarloSmallSampleSize;
                end
                geneIsContained = zeros(size(randGenePoints,1),1);
                montePointIsInside = zeros(monteCarloSampleSize,1);

                for region = 1 : numel(numVertices)
                    confidenceCont = zeros(numVertices(region),2);
                    for r = 1 : numVertices(region)
                        confidenceCont(r,1) = densityContPoints(1,confidenceContInd(region)+r);
                        confidenceCont(r,2) = densityContPoints(2,confidenceContInd(region)+r);
                    end 
                    % estimate contour:grid area ratio
                    montePointIsInside = montePointIsInside +...
                        inpolygon(randPoints(:,1),randPoints(:,2),...
                        confidenceCont(:,1),confidenceCont(:,2));

                    % convert contour (Grid) points back to original scale (TPM)
                    confidenceCont(:,1) = (confidenceCont(:,1) - 1)*gridwidth + startx - gridwidth;
                    confidenceCont(:,2) = (confidenceCont(:,2) - 1)*gridheight + starty - gridheight;

                    % find % of gene prob. contained by contour
                    geneIsContained  = geneIsContained +...
                        inpolygon(randGenePoints(:,1),randGenePoints(:,2),...
                        confidenceCont(:,1),confidenceCont(:,2));
                end

                contGridAreaRatio = numel(find(montePointIsInside))/...
                    monteCarloSampleSize;
                geneProbContained = numel(find(geneIsContained))/...
                    size(geneIsContained,1);


    %             if contGridAreaRatioCutoff > contGridAreaRatio &&...
    %                     geneProbContained > geneProbContainedCutoff
    %                 contourIsNotCorrect = 0;
    %             end

                if testingContourSelection
                    geneProbContainedInLv(lv) = geneProbContained;
                    areaRatioOfLv(lv) = contGridAreaRatio;
                    if lv >= numel(lvls)
                        testingContourSelection = 0;
                    end
                end
            else
                zeroVertices = 1;
            end


            % return non working pair
            if contourIsNotCorrect && lv >= numel(lvls)
                contourIsNotCorrect = 0;
                zeroVertices = 0;
            end
        end

        % find optimal contour indices
        contEffectiveness = (1 - abs(targetContainment-geneProbContainedInLv))...
            ./log2(areaRatioOfLv+2);

        % find most effective contour
        [~,optimalContLv] = max(contEffectiveness);
        validateLv = optimalContLv;
        trueOptimalContLv = optimalContLv;
        if validation
            if validateAllContours
                validateLvls = 1:numel(lvls);
            else
                validateLvls = sort(unique(...
                    [optimalContLv-2,optimalContLv-1,optimalContLv,optimalContLv+1,optimalContLv+2,...
                    1:round(numel(lvls)/validateN):numel(lvls)]));
            end
            validateOS = zeros(numGene,numel(validateLvls));
            %validateOgenes = zeros(numel(validateLvls),1);
            validateLv = numel(validateLvls);
            validateRefContour = find(validateLvls == optimalContLv);
        end
        
        
        % save optimal contour info (containment/size/effectivness)
        contourLoopsOptimalInfo(c,1) = geneProbContainedInLv(optimalContLv);
        contourLoopsOptimalInfo(c,2) = areaRatioOfLv(optimalContLv);
        contourLoopsOptimalInfo(c,3) = contEffectiveness(optimalContLv);
        if c > 1
            if contourLoopsOptimalInfo(c,3) <= contourLoopsOptimalInfo(c-1,3)
                contourLoopsOptimalInfo = contourLoopsOptimalInfo(1:c,:);
                c = contourLoopsMax;
            end
        end
        
        if outputPlots && c == contourLoopsMax
            hold off
            axis equal
            title('Contours');
            xticks([]);
            yticks([]);

            ax4 = nexttile;
            ax5 = nexttile;
            ax6 = nexttile;
            if validation
                clr = parula(size(validateOS,2));
                clr(validateRefContour,:) = [1 0 0];
                figure;
                axValContours = gca;
                hold(axValContours,'on')
            end
        end
        
        % select contour height range
        if optimalContLv <= countourCutoff
            contourRange = [lvls(1):...
                abs(lvls(optimalContLv + countourCutoff)-lvls(1))/numContours:...
                lvls(optimalContLv + 1)];
        elseif optimalContLv >= numel(lvls) - countourCutoff
            contourRange = [lvls(optimalContLv - countourCutoff):...
                abs(lvls(optimalContLv)-lvls(optimalContLv - countourCutoff))/numContours:...
                lvls(optimalContLv)];
        else
            contourRange = [lvls(optimalContLv - countourCutoff):...
                abs(lvls(optimalContLv + countourCutoff)-...
                lvls(optimalContLv - countourCutoff))/numContours:...
                lvls(optimalContLv + countourCutoff)];
        end
    end
    
    disp('Assigning outlier scores')
    
    while validateLv > 0
        if validation
            optimalContLv = validateLvls(validateLv);
            disp(['validating contour level: ',num2str(validateLvls(validateLv))])
        end
        confidenceContInd = zeros(size(densityContPoints,2),1);
        temp = 1;
        while temp < size(densityContPoints,2)
            if densityContPoints(1,temp) == lvls(optimalContLv)
                confidenceContInd(temp) = 1;
            end
            temp = temp + densityContPoints(2,temp) + 1;
        end
        confidenceContInd = find(confidenceContInd);
        numVertices = densityContPoints(2,confidenceContInd);

        
        OutlierScore = zeros(numGene,1);
        OutlierDistScore = zeros(numGene,1);
        firstPointCalc = ones(numGene,1);
        printConfidenceCont = zeros(1,2);
        randGene = randi(numGene);
        monteCarloSampleSize = 1000;

        for region = 1 : numel(numVertices)
            confidenceCont = zeros(numVertices(region),2);
            for r = 1 : numVertices(region)
                confidenceCont(r,1) = densityContPoints(1,confidenceContInd(region)+r);
                confidenceCont(r,2) = densityContPoints(2,confidenceContInd(region)+r);
            end
            % convert contour (Grid) points back to original scale (TPM)
            confidenceCont(:,1) = (confidenceCont(:,1) - 1)*gridwidth + startx - gridwidth;
            confidenceCont(:,2) = (confidenceCont(:,2) - 1)*gridheight + starty - gridheight;

            % save points from every region separated by [0 0]
            printConfidenceCont = [printConfidenceCont;confidenceCont;[0 0]];

            % find area of gene probability outside of
            % confidence contour (monte carlo area est.)
            for k = 1 : numGene

                % random points within 2sigma
                randPoints = [geneMeanX(k) + geneSTDX(k)*randn(monteCarloSampleSize,1),...
                    geneMeanY(k) + geneSTDY(k)*randn(monteCarloSampleSize,1)];
                OutlierScore(k) = OutlierScore(k) +...
                    numel(find(inpolygon(randPoints(:,1),randPoints(:,2),...
                    confidenceCont(:,1),confidenceCont(:,2))))/monteCarloSampleSize;
                
                % find shortest distance to contour of mean point
                for contourSeg = 1 : size(confidenceCont,1) - 1
                    %initialize w/ z-direction for cross function
                    tempA = zeros(3,1);
                    tempB = zeros(3,1);
                    
                    tempA(1:2) = confidenceCont(contourSeg,:) -...
                        confidenceCont(contourSeg+1,:);
                    tempB(1:2) = [geneMeanX(k),geneMeanY(k)] -...
                        confidenceCont(contourSeg+1,:);
                    
                    distToSeg = norm(cross(tempA,tempB))/norm(tempA);
                    
                    if distToSeg < abs(OutlierDistScore(k)) || firstPointCalc(k)
                        firstPointCalc(k) = 0;
                        if inpolygon(geneMeanX(k),geneMeanY(k),...
                                confidenceCont(:,1),confidenceCont(:,2))
                            OutlierDistScore(k) = -1*distToSeg;
                        else
                            OutlierDistScore(k) = distToSeg;
                        end
                    end
                end
                if outputPlots && k == randGene
                    scatter(ax4,randPoints(:,1),randPoints(:,2),25,'r.');
                    hold(ax4,'on')
                    title(ax4,'Scoring Single Gene');
                    xlabel(ax4,'TPM');
                    ylabel(ax4,'TPM');
                    scatter(ax4,geneMeanX(k),geneMeanY(k),50,'bx','LineWidth',2);
                    plot(ax4,confidenceCont(:,1),confidenceCont(:,2),'k');
                    legend(ax4,{'MC points','mean',''})
                    
                    if validation &&...
                            (abs(validateLv - validateRefContour) > 2 ||...
                            validateLv == validateRefContour)
                        plot(axValContours,confidenceCont(:,1),...
                            confidenceCont(:,2),'LineWidth',2,...
                        'Color',clr(validateLv,:));
                    end
                end
                %[region k]
            end
        end
        
        if removeHighLowExpr
            % remove low/high expression outlier Scores
            meanDataXSort = sort(dataX);
            meanDataYSort = sort(dataY);
            lowX = meanDataXSort(floor(...
                (1-contourLoopsOptimalInfo(end,1))/2*numGene));
            lowY = meanDataYSort(floor(...
                (1-contourLoopsOptimalInfo(end,1))/2*numGene));
            highX = meanDataXSort(ceil(...
                (contourLoopsOptimalInfo(end,1)+...
                (1-contourLoopsOptimalInfo(end,1))/2)*numGene));
            highY = meanDataYSort(ceil(...
                (contourLoopsOptimalInfo(end,1)+...
                (1-contourLoopsOptimalInfo(end,1))/2)*numGene));
        
        
            for i = 1 : numGene
                if (mean(dataX(i,:)) <= lowX && mean(dataY(i,:)) <= lowY) ||...
                        (mean(dataX(i,:)) >= highX && mean(dataY(i,:)) >= highY)
                    OutlierScore(i) = 1;
                end
            end
        end
        
        % set outlier score max to 1
        OutlierScore(OutlierScore>1) = 1;
        
        if validation
            validateOS(:,validateLv) = OutlierScore;
            validateLv = validateLv - 1;
        else
            validateLv = 0;
        end
    end
    
    if validation
        OutlierScore = validateOS(:,validateRefContour);
        optimalContLv = trueOptimalContLv;
    end
    
    % rank order distScores < 0
    Dpenalty = zeros(numel(OutlierDistScore),1);
    %innerDist = OutlierDistScore(OutlierDistScore < 0);
    [~,Dpenalty(OutlierDistScore < 0)] = sort(abs(OutlierDistScore(OutlierDistScore < 0)));
    [~,Dpenalty(OutlierDistScore < 0)] = sort(Dpenalty(OutlierDistScore < 0));
    Dpenalty = Dpenalty/numel(find(OutlierDistScore < 0));
    
    % ---adjust OS by STD---
    OutlierScore = 1 - OutlierScore - Dpenalty;
    OutlierScore(OutlierScore < 0) = 0;
    
    % rank order STD
%     [~,geneSTDXrank] = sort(geneSTDX);
%     [~,geneSTDYrank] = sort(geneSTDY);
%     [~,geneSTDXrank] = sort(geneSTDXrank);
%     [~,geneSTDYrank] = sort(geneSTDYrank);
%     for k = 1 : numGene
%         OutlierScore(k) = OutlierScore(k)*...
%             
%             %(1-min([totalSTD,1])*sqrt((geneSTDXrank(k)^2+geneSTDYrank(k)^2)/2)/numGene);
%             %(1-(geneSTDXrank(k)+geneSTDYrank(k))/(2*numGene));
%             %(sqrt((geneSTDX(k)+geneSTDY(k))/2));
%     end
            
    % display outliers w/ classification above 0.95 outlier score
    [~,inCount] = sort(OutlierScore);
    in = false(size(dataY,1),1);
    in(inCount(1:floor(0.95*size(dataY,1)))) = 1;
    in = find(in);

    out = [1:size(dataY,1)]';
    out = out(~ismember(out,in));
            

    if outputPlots
        disp('Displaying results')
        hold(ax4,'off')
    end

    % display plots
    if outputPlots
        scatter(ax5,geneMeanX(in),geneMeanY(in),25,'c.');
        hold(ax5,'on')
        scatter(ax5,geneMeanX(out),geneMeanY(out),25,'r.');
        stopPoints = intersect(find(printConfidenceCont(:,1) == 0),...
            find(printConfidenceCont(:,2) == 0));
        for v = 1 : numel(stopPoints) - 1
            plot(ax5,printConfidenceCont(stopPoints(v)+1:stopPoints(v+1)-1,1),printConfidenceCont(stopPoints(v)+1:stopPoints(v+1)-1,2),'g','LineWidth',3);
        end
        hold(ax5,'off')
        axis equal
        title(ax5,'Top 5% OS Genes');
        xlabel(ax5,'TPM');
        ylabel(ax5,'TPM');

        histogram(ax6,OutlierScore,'FaceColor','b');
        set(gca, 'YScale', 'log')
        title('OS Distribution');
        ylabel('Gene Count');
        xlabel('OS');
        
        figure;
        surf(log2(densityMat/numGene+1),'FaceColor','interp',...
            'FaceLighting','gouraud','EdgeColor','none')
        title('Expression Density Topology')
        xlabel('Sample 1 TPM')
        ylabel('Sample 2 TPM')
        zlabel('Cumulative PDF')
        xticklabels(round((xticks - 1)*gridwidth + startx - gridwidth,1))
        yticklabels(round((yticks - 1)*gridwidth + startx - gridwidth,1))
        
        if validation
            validateOS = 1 - validateOS;
            
            % update contour plot
            hold(axValContours,'off')
            xticks(axValContours,[])
            yticks(axValContours,[])
            xlabel(axValContours,[])
            ylabel(axValContours,[])
            title(axValContours,'Selected Contours')
            set(axValContours,'LineWidth',2,'FontWeight','bold')
            
            % find rank of OS in each contour level
            clustOS = validateOS(randperm(numGene,min([numGene,1000])),:);
            [~,clustOSRank] = sort(clustOS,'ascend');
            [~,clustOSRank] = sort(clustOSRank);
            [~,refInd] = sort(clustOSRank(:,validateRefContour));
            
%             clustergram(clustOS,'Symmetric','false','Colormap','turbo',...
%                 'Cluster','column');
            
            figure;
            imagesc(clustOSRank(refInd,:));
            title('OS Rank Sorted By Optimal Contour')
            xticks(1:numel(validateLvls))
            xticklabels(validateLvls)
            xlabel('Contour Level')
            ylabel('Rank Ordered Genes')
            ax = gca;
            set(ax,'LineWidth',2,'FontWeight','bold')
            
            figure;
            ax = gca;
            hold(ax,'on')
            figure;
            ax2 = gca;
            hold(ax2,'on')
            for i = 1 : size(validateOS,2)
                if i <= validateRefContour
                    plot(ax,sort(validateOS(:,i)),'LineWidth',3,...
                        'Color',clr(i,:));
                end
                if i >= validateRefContour
                    plot(ax2,sort(validateOS(:,i)),'LineWidth',3,...
                        'Color',clr(i,:));
                end
            end
            plot(ax2,sort(validateOS(:,validateRefContour)),...
                'LineWidth',3,'Color',clr(validateRefContour,:));
            hold(ax,'off')
            lgnd = string(validateLvls);
            lgnd(validateRefContour) = "O.C.";
            lgd = legend(ax,lgnd(1:validateRefContour),'Location','eastoutside');
            title(lgd,'contour level')
            lgd = legend(ax2,lgnd(validateRefContour:end),'Location','eastoutside');
            title(lgd,'contour level')
            xlabel(ax,'Rank ordered genes')
            ylabel(ax,'Outlier Score')
            xlabel(ax2,'Rank ordered genes')
            ylabel(ax2,'Outlier Score')
            set(ax,'LineWidth',2,'FontWeight','bold')
            set(ax2,'LineWidth',2,'FontWeight','bold')
            
            [~,validateOSsortInd] = sort(validateOS,'descend');
            
            % OS difference between highest 0-1% genes and 4-5% genes
            deltaPercents = [0,1; 4,5];
            figure;
            OSdelta = zeros(size(validateOS,2),1);
            deltaInd = round((deltaPercents/100)*numGene);
            deltaInd(deltaInd == 0) = 1;
            for i = 1 : numel(OSdelta)
                OSdelta(i) = mean(validateOS(validateOSsortInd(deltaInd(1,1):deltaInd(1,2),i),i)) -...
                    mean(validateOS(validateOSsortInd(deltaInd(2,1):deltaInd(2,2),i),i));
            end
            %plot(OSdelta,'LineWidth',3)
            b = bar(OSdelta);
            b.FaceColor = 'flat';
            b.CData = clr;
            xticks(1:numel(OSdelta))
            xticklabels(validateLvls)
            xlabel('Contour level')
            ylabel("\Delta OS ("+num2str(deltaPercents(1,1))+" - "+num2str(deltaPercents(1,2))+...
                "% : "+num2str(deltaPercents(2,1))+" - "+num2str(deltaPercents(2,2))+"%)")
            set(gca,'LineWidth',2,'FontWeight','bold')
            
            
            numOutlier = zeros(size(validateOS,2),1);
            outlierOverlap = zeros(size(validateOS,2),1);
            outlierTopOverlap = zeros(size(validateOS,2),1);
            for i = 1 : numel(numOutlier)
                numOutlier(i) = numel(find(validateOS(:,i) >= validationOScutoff));
                outlierOverlap(i) = numel(intersect(...
                    find(validateOS(:,i) >= validationOScutoff),...
                    find(validateOS(:,validateRefContour) >= validationOScutoff)));
                outlierTopOverlap(i) = numel(intersect(...
                    validateOSsortInd(1:round(0.05*numGene),i),...
                    validateOSsortInd(1:round(0.05*numGene),validateRefContour)));
            end
            outlierOverlap = outlierOverlap/numOutlier(validateRefContour);
            outlierTopOverlap = outlierTopOverlap/round(0.05*numGene);
            
            figure;
            subplot(4,1,1)
            b = bar(numOutlier/numGene);
            b.FaceColor = 'flat';
            b.CData = clr;
            %b.CData(validateRefContour,:) = [1 0 0];
            text(validateRefContour - 0.5,numOutlier(validateRefContour)+0.2*max(numOutlier),...
                {'Optimal','Contour','\downarrow'},'Color','red')
            xticklabels(validateLvls)
            xlabel('Contour Level')
            ylabel('TAG %')
            title(['Genes w/ OS > ',num2str(validationOScutoff)])
            ax = gca;
            set(ax,'LineWidth',2,'FontWeight','bold')
            
            subplot(4,1,2)
            b = bar(outlierOverlap);
            b.FaceColor = 'flat';
            %b.CData(validateRefContour,:) = [1 0 0];
            b.CData = clr;
            text(validateRefContour - 0.5,numOutlier(validateRefContour)+0.2*max(numOutlier),...
                {'Optimal','Contour','\downarrow'},'Color','red')
            xticklabels(validateLvls)
            xlabel('Contour Level')
            ylabel('% Agreement')
            title(['TAG Agreement w/ Optimal Contour: OS >',...
                num2str(validationOScutoff)])
            ax = gca;
            set(ax,'LineWidth',2,'FontWeight','bold')
            
            subplot(4,1,3)
            b = bar(outlierTopOverlap);
            b.FaceColor = 'flat';
            %b.CData(validateRefContour,:) = [1 0 0];
            b.CData = clr;
            text(validateRefContour - 0.5,numOutlier(validateRefContour)+0.2*max(numOutlier),...
                {'Optimal','Contour','\downarrow'},'Color','red')
            xticklabels(validateLvls)
            xlabel('Contour Level')
            ylabel('% Agreement')
            title('Top 5% OS Agreement w/ Optimal Contour')
            ax = gca;
            set(ax,'LineWidth',2,'FontWeight','bold')
            
            subplot(4,1,4)
            boxplot(validateOS);
            xticklabels(validateLvls)
            xlabel('Contour Level')
            ylabel('OS')
            title('OS Scores per Contour Level')
            ax = gca;
            set(ax,'LineWidth',2,'FontWeight','bold')
        end
    
    
        figure;
        subplot(3,1,1)
        plot(geneProbContainedInLv,'LineWidth',3)
        yline(targetContainment,'--r','T: target');
        hold on
        plot(optimalContLv,geneProbContainedInLv(optimalContLv),'r*','LineWidth',2,'MarkerSize',10)
        text(optimalContLv+2,geneProbContainedInLv(optimalContLv),'Optimal Contour','Color','red')
        hold off
        title('Contour Probability Coverage')
        xlabel('Contour Level')
        ylabel('% of Prob')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        subplot(3,1,2)
        plot(areaRatioOfLv,'LineWidth',3)
        hold on
        plot(optimalContLv,areaRatioOfLv(optimalContLv),'r*','LineWidth',2,'MarkerSize',10)
        text(optimalContLv+2,areaRatioOfLv(optimalContLv),'Optimal Contour','Color','red')
        hold off
        title('Contour Size')
        xlabel('Contour Level')
        ylabel('log_2(S + 1)')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        subplot(3,1,3)
        plot(contEffectiveness,'LineWidth',3)
        hold on
        plot(optimalContLv,contEffectiveness(optimalContLv),'r*','LineWidth',2,'MarkerSize',10)
        text(optimalContLv+2,contEffectiveness(optimalContLv),'Optimal Contour','Color','red')
        hold off
        title('Contour Effectiveness')
        xlabel('Contour Level')
        ylabel('CE')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        sgtitle('Optimal Contour Selection Criteria')
        
        figure;
        subplot(3,1,1)
        plot(contourLoopsOptimalInfo(:,1),'LineWidth',3)
        yline(targetContainment,'--r','T: target');
        title('Contour Probability Coverage')
        xlabel('Iteration')
        ylabel('% of Prob')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        subplot(3,1,2)
        plot(contourLoopsOptimalInfo(:,2),'LineWidth',3)
        title('Contour Size')
        xlabel('Iteration')
        ylabel('log_2(S + 1)')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        subplot(3,1,3)
        plot(contourLoopsOptimalInfo(:,3),'LineWidth',3)
        title('Contour Effectiveness')
        xlabel('Iteration')
        ylabel('CE')
        ax = gca;
        set(ax,'LineWidth',2,'FontWeight','bold')
        sgtitle('Contour Optimization')
        figure; % catch new contour plots
    end
end

