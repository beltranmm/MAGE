classdef MAGE_functions
    
    methods (Static)
        
        %% simData
%         DESCRIPTION: Generate simulated data for testing of MAGE.
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
%               sim_profile -- double -- (gene, samples)
%
        
        function sim_profile = simData(m,n,mu,sigma,DErate)
            
            sim_profile = zeros(m,2*n);
            for i = 1 : m
                geneMean = abs(mu + mu*randn(1));
                geneSD = abs(sigma + sigma*randn(1));
                
                for j = 1 : n
                    sim_profile(i,j) = geneMean + geneSD*randn(1);
                end
                
                if DErate >= rand
                    geneMean = abs(mu + mu*randn(1));
                end
                geneSD = abs(sigma + sigma*randn(1));
                
                for j = 1 : n
                    sim_profile(i,j+n) = geneMean + geneSD*randn(1);
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

    end
end