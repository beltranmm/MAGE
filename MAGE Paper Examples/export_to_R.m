%% export data to R (.csv) (11/1/2023)
% After profile has been loaded workspace should contain:
% profile ---> contains un-normalized RNA-seq read counts
% geneName ---> list of gene names
% sampleName ---> list of sample names

if size(sampleName,1) > 0
    sampleName = sampleName';
end

controlBatch = 1;
treatmentBatch = 1;
batch = zeros(numel(sampleName),1);
condition = zeros(numel(sampleName),1);

for i = 1 : numel(sampleName)
    if strcmp(sampleName(i),'control')
        condition(i) = 1;
        batch(i) = controlBatch;
        controlBatch = controlBatch + 1;
    else
        condition(i) = 2;
        batch(i) = treatmentBatch;
        treatmentBatch = treatmentBatch + 1;
    end
end

writematrix(profile,strcat(dataDir,'edgeR\cts.csv'));
writematrix(sampleName,strcat(dataDir,'edgeR\colData.csv'));
writematrix(geneName,strcat(dataDir,'edgeR\geneName.csv'));
writematrix(batch,strcat(dataDir,'edgeR\batch.csv'));
writematrix(condition,strcat(dataDir,'edgeR\condition.csv'));