% 


function [OutlierScore, FDR] = MAGE(dataX,dataY,varargin)    
    
    %% check input parameters
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
    
    
    
    %% calculate gene statistics
    numGene = size(dataY,1);
    
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
    
    
    
    %% fit grid
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
    
    
    
    %% sum gene probabilities
    disp('Summing gene probabilities')
    
    % find density within each grid region
    for k = 1 : numGene
        % start at 1,1 in grid
        gridLocationX = startx;
        for gridIndX = 1 : gridDensity

            gridLocationY = starty;
            for gridIndY = 1 : gridDensity
                
                % determine densityMat in nonnegative region
                if gridLocationX > 0 && gridLocationY > 0 
                    densityMat(gridIndY,gridIndX) =...
                        densityMat(gridIndY,gridIndX) +...
                        (normpdf(gridLocationX,geneMeanX(k),geneSTDX(k)) *...
                        normpdf(gridLocationY,geneMeanY(k),geneSTDY(k)));
                end

                gridLocationY = gridLocationY + gridheight;
            end
            gridLocationX = gridLocationX + gridwidth;
        end
    end


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

    
    
    %% select CER
    disp('Determining characteristic expression region ')
    
    % pad boundary to ensure closed contour lines
    densityMat = padarray(densityMat,[1 1], 0, 'both');
    
    % set negative space to zero
    
    
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

                % Monte-Carlo sampling from CPDF
                monteCarloSampleSize = 1000;
                randPoints = [gridDensity*rand(monteCarloSampleSize,1),...
                    gridDensity*rand(monteCarloSampleSize,1)];
                monteCarloSmallSampleSize = 10;
                randGenePoints = zeros(numGene*monteCarloSmallSampleSize,2);
                temp = 1;
                for k = 1 : numGene
                    for smallSample = 1 : monteCarloSmallSampleSize
                        randGenePoints(temp,:) =...
                            [geneMeanX(k) + geneSTDX(k)*randn(1,1),...
                            geneMeanY(k) + geneSTDY(k)*randn(1,1)];
                        temp = temp + 1;
                    end
                end
                
                % keep only points from 1st quadrant
                randGenePoints = randGenePoints(randGenePoints(:,1) > 0 & randGenePoints(:,2) > 0,:);
                
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

                geneProbContained = numel(find(geneIsContained))/...
                    size(geneIsContained,1);

                if testingContourSelection
                    geneProbContainedInLv(lv) = geneProbContained;
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
        contEffectiveness = 1 - abs(targetContainment-geneProbContainedInLv);

        % find most effective contour
        [~,optimalContLv] = max(contEffectiveness);
        
        % save optimal contour info (containment/size/effectivness)
        contourLoopsOptimalInfo(c,1) = geneProbContainedInLv(optimalContLv);
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
    
    
    
    %% determine OS_raw
    disp('Assigning outlier scores')
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
            
            % only count points in 1st quadrant
            randPoints = randPoints(randPoints(:,1) > 0 & randPoints(:,2) > 0,:);
            
            OutlierScore(k) = OutlierScore(k) +...
                numel(find(inpolygon(randPoints(:,1),randPoints(:,2),...
                confidenceCont(:,1),confidenceCont(:,2))))/size(randPoints,1);

            % find shortest distance to contour of mean point
            for contourSeg = 1 : size(confidenceCont,1) - 1
                
                tempA = confidenceCont(contourSeg,:) -...
                    confidenceCont(contourSeg+1,:);
                tempB = [geneMeanX(k),geneMeanY(k)] -...
                    confidenceCont(contourSeg+1,:);

                distToSeg = norm(tempB - tempA);

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
            end
            %[region k]
        end
    end

    
    
    %% OS distance adjustment
    if removeHighLowExpr
        % remove low/high expression outlier Scores
        meanDataXSort = sort(mean(dataX,2));
        meanDataYSort = sort(mean(dataY,2));
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
    
    
    % rank order distScores < 0
    Dpenalty = zeros(numel(OutlierDistScore),1);
    [~,Dpenalty(OutlierDistScore < 0)] = sort(abs(OutlierDistScore(OutlierDistScore < 0)));
    [~,Dpenalty(OutlierDistScore < 0)] = sort(Dpenalty(OutlierDistScore < 0));
    Dpenalty = Dpenalty/numel(find(OutlierDistScore < 0));
    
    % ---adjust OS by STD---
    OutlierScore = 1 - OutlierScore - Dpenalty;
    OutlierScore(OutlierScore < 0) = 0;
    
    
    
    %% display        
    % display genes with highest 0.05 of OS values
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
        scatter(ax5,geneMeanX(in),geneMeanY(in),75,[0.7,0.7,0.7],'.');
        hold(ax5,'on')
        scatter(ax5,geneMeanX(out),geneMeanY(out),75,'r.');
        stopPoints = intersect(find(printConfidenceCont(:,1) == 0),...
            find(printConfidenceCont(:,2) == 0));
        for v = 1 : numel(stopPoints) - 1
            plot(ax5,printConfidenceCont(stopPoints(v)+1:stopPoints(v+1)-1,1),printConfidenceCont(stopPoints(v)+1:stopPoints(v+1)-1,2),'k','LineWidth',3);
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
    
        figure; % catch new contour plots
    end
    
    
    %% Calculate FDR
    if nargout >= 2
        
        disp('Calculating FDR')
        
        % permuate samples
        dataTotal = [dataX, dataY];
        for i = 1 : 10
            dataTotal = dataTotal(:,randperm(size(dataTotal,2))');
        end


        % Keep only equally expressed (EE) genes
        % by removing genes with highest 5% OS
        [~,OS_sortInd] = sort(OutlierScore);
        dataTotal = dataTotal(OS_sortInd(1:round(0.95*numel(OS_sortInd))),:);
        dataX = dataTotal(:,1:size(dataX,2));
        dataY = dataTotal(:,1:size(dataY,2));
        
        % recursive call with permutated data
        OutlierScore_perm = MAGE(dataX,dataY,gridDensity,...
            numContours,outputPlots,targetContainment,...
            removeHighLowExpr,contourLoopsMax,numStartingContours);
        
        %FDR adj.
        OutlierScore_perm(OutlierScore_perm > 0.8) = 0;
        
        % calculate FDR at each OS
        OS_FDR = zeros(100,2);
        outlierScoreCutoff = 0;
        
        for i = 1 : size(OS_FDR,1)
            numAbove = numel(find(OutlierScore >= outlierScoreCutoff));
            numAbove_P = numel(find(OutlierScore_perm >= outlierScoreCutoff));
            
            OS_FDR(i,1) = outlierScoreCutoff;
            OS_FDR(i,2) = numAbove_P/(numAbove + numAbove_P);
            
            outlierScoreCutoff = outlierScoreCutoff + 1/size(OS_FDR,1);
        end
        OS_FDR(isnan(OS_FDR(:,2)),2) = 0;
        
        figure;
        scatter(OS_FDR(:,1),OS_FDR(:,2),'filled')
        xlabel('OS')
        ylabel('FDR')
        set(gca,'TickDir','out');
        
        
        
        % assign FDR to genes
        FDR = zeros(size(OutlierScore,1),1);
        for i = 1 : numel(FDR)
            [~,tmp] = min(abs(OutlierScore(i) - OS_FDR(:,1)));
            FDR(i) = OS_FDR(tmp,2);
        end 
    end
    
    
end

