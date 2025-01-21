% Transfer from workspace to csv
% Matthew Beltran
% 1/20/2025

% gt3
%fileName = 'gt3_NCBI.csv';
%profileVar = profile_GDS4059fullabemaciclib;
%geneNameVar = geneName_GDS4059fullabemaciclib;
%sampleNameVar = sampleName_GDS4059fullabemaciclib;

% mtor
fileName = 'mtor_NCBI.csv';
profileVar = profile;
geneNameVar = geneName;
sampleNameVar = sampleName;

profile = cell(numel(geneNameVar)+1,numel(sampleNameVar)+1);
for i = 1 : numel(sampleNameVar)
    profile(1,i+1) = cellstr(sampleNameVar(i));
end

for i = 1 : numel(geneNameVar)
    profile(1+i,1) = cellstr(geneNameVar(i));
end

for i = 1 : numel(geneNameVar)
    for j = 1 : numel(sampleNameVar)
        profile(1+i,j+1) = num2cell(profileVar(i,j));
    end
end


profile = table(profile);

writetable(profile,fileName);

