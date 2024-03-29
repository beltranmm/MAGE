classdef MAGE_functions
    
    methods (Static)
        
        %% simData
%         DESCRIPTION: Generate simulated data for testing of MAGE.
%         
%         INPUT:
%         
%               m -- double -- number of genes
%
%               n -- double -- number of sample replicates
%
%               mu -- double -- mean expression of all genes
%
%               sigma  -- double -- mean standard deviation of all genes
%
%               DErate -- double -- rate of differential expression
%
%               
%         OUTPUT:
%
%               sim_profile -- double -- (gene, samples)
%               
%               DEind -- double -- indices of differntially expressed genes

        
        function [sim_profile,DEind] = simData(m,n,mu,sigma,DErate)
            
            sim_profile = zeros(m,2*n);
            DEind = false(m,1);
            for i = 1 : m
                geneMean = abs(mu + mu*randn(1));
                geneSD = abs(sigma + sigma*randn(1));
                
                for j = 1 : n
                    sim_profile(i,j) = abs(geneMean + geneSD*randn(1));
                end
                
                if DErate >= rand
                    geneMean = abs(mu + mu*randn(1));
                    DEind(i) = true;
                end
                geneSD = abs(sigma + sigma*randn(1));
                
                for j = 1 : n
                    sim_profile(i,j+n) = abs(geneMean + geneSD*randn(1));
                end
            end
            
        end
        
%% addNoise
%         DESCRIPTION: Add normally distributed random noise to 
%         expression profile.
%         
%         INPUT:
%         
%               profile_in -- double -- number of genes
%
%               shift -- double -- mean expression of all genes
%
%               sigma  -- double -- standard deviation of normally
%                                   distributed noise
%               
%         OUTPUT:
%
%               profile_out -- double -- (gene, samples)
%               

        
        function profile_out = addNoise(profile_in,shift,sigma)
            
            profile_out = zeros(size(profile_in));
            
            for i = 1 : size(profile_in,1)
                for j = 1 : size(profile_in,2)
                    profile_out(i,j) = profile_in(i,j) +...
                        sigma(j)*(shift(j) + randn(1));
                end
            end
            
        end
        
        %% fltr
%         DESCRIPTION: Filter out genes with low expression. If a gene is
%         expressed below "level" in more than "sampleMin" number of
%         samples then it is removed from the filtered profile
%         
%         INPUT:
%         
%               profile -- double -- (genes, samples)
%
%               level -- double -- minimum expression level
%
%               sampleMin -- double -- minimum number of samples needed to
%            	be above minimum level
%         
%         OUTPUT:
%
%               filtered_profile -- double -- (gene, samples)
%               
%               filteredInd -- double -- The index from the original
%               profile for every gene in the filyered profile
%
        
        function [filtered_profile, filteredInd] = fltr(profile,level,sampleMin) 
            keepGene = false(size(profile,1),1);
            for g = 1 : size(profile,1)
                keepGene(g) =...
                    numel(find(profile(g,:) > level)) >= sampleMin;
            end
            filtered_profile = profile(keepGene,:);
            filteredInd = keepGene;
        end
        
        %% ROC
%         DESCRIPTION: Calculate the true/false positive rates based on
%         given test statistic by sweeping classification cutoff values
%         
%         INPUT:
%         
%               t -- double -- test statistic values
%
%               trueInd -- double -- indices of true target classifications
%         
%         OUTPUT:
%
%               TPR -- double -- points for true-positive rate
%
%               FPR -- double -- points for false-positive rate
%
        function [TPR,FPR] = ROC(t,trueInd)
            numpts = 1000;
            
            step = (max(t) - min(t))/numpts;
            cutoff = min(t);
            
            TPR = zeros(numpts,1);
            FPR = zeros(numpts,1);

            for i = 1 : numel(TPR)
                TPR(i) = numel(find(t(trueInd) >= cutoff))/numel(find(trueInd));
                FPR(i) = 1 - numel(find(t(~trueInd) < cutoff))/numel(find(~trueInd));

                cutoff = cutoff + step;
            end
        end
        
        %% AUC
%         DESCRIPTION: Calculate area under the curve from ROC using
%         monte-carlo area est.
%         
%         INPUT:
%         
%               FPR -- double -- X-axis points from ROC (false-positives)
%
%               TPR -- double -- Y-axis points from ROC (true-positives)
%         
%         OUTPUT:
%
%               auc -- double -- area under the curve
%
        function auc = AUC(FPR,TPR)
            numMCpts = 1000000;
            
            % enclose area outlined by ROC
            c = [FPR,TPR;[0,0];[1,0]];
            c = unique(c,'stable','rows');
            
            rp = rand(numMCpts,2);
            auc = numel(find(inpolygon(rp(:,1),rp(:,2),c(:,1),c(:,2))))/numMCpts;
        end

    end
end