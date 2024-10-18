%% setup

% select data directory
dataDir = strcat(pwd,'\workspaces\');

% load contour functions
f = DEG_contour_functions;


profile_complete = zeros();

%% DBSCAN: breast cancer gamma T3 treatment/control
load(strcat(dataDir,'breast_cancer_abemaciclib_NCBI.mat'));

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

title1 = "breast cancer (abemaciclib)";

%f.filterAnalysis(profile,[1,3,5,8],[1,2,3,6]);
[~, ind] = f.fltr(profile,1,6);
profile = profile(ind,:);
geneName = geneName(ind);

% log shift
%profile = f.logNorm(profile,1,0);

profile_complete = profile;
geneName_complete = geneName;

%------------run DBSCAN--------------------------

n = size(profile,1);

% best parameters for full g-T3 profile
r = 0.3;
minPts = 200;
OScutoff = 0.65;
excludeTop = 0.05;
title_tmp = "\gamma - T3";
xL_tmp = "Untreated MCF-7 cells (log_2(TPM+1))";
yL_tmp = "\gamma-T3 MCF-7 cells (log_2(TPM+1))";

runTopo = 1;
gridDensity = 100;
numContours = 20;
%--------------------------------------------------------------------------
sampleClass = ["control";"treatment"];



profileReduced = [mean(profile(:,strcmp(sampleName,sampleClass(1))),2),...
    mean(profile(:,strcmp(sampleName,sampleClass(2))),2)];

nOut = zeros(numel(r),numel(minPts));
nClusts = zeros(numel(r),numel(minPts));

for i = 1 : numel(r)
    for j = 1 : numel(minPts)
        dbOutput = dbscan(profileReduced(1:n,:),r(i),minPts(j));
        nOut(i,j) = sum(dbOutput == -1);
        nClusts(i,j) = numel(unique(dbOutput)) - 1;
        [i j]
    end
end


grp = "clust "+string(dbOutput);
grp(dbOutput == -1) = "outlier";

figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,[],[],10)
title("DBSCAN clusters")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

grp(dbOutput ~= -1) = "clustered";

figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,'cr',[],10)
title("DBSCAN classification")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

profileReduced_sorted = sort(profileReduced,1,'descend');
excludeCutoffs = profileReduced_sorted(round(n*excludeTop),:);
grp(profileReduced(:,1) > excludeCutoffs(1) &...
    profileReduced(:,2) > excludeCutoffs(2)) = "clustered";
%grp(profileReduced(:,2) > excludeCutoffs(2)) = "clustered";
figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,'cr',[],10)
xline(excludeCutoffs(1),'m--','LineWidth',2)
yline(excludeCutoffs(2),'m--','LineWidth',2)
title("DBSCAN classification (exluding top "+excludeTop*100+"%)")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
legend({"clustered","outlier"}) %#ok<CLARRSTR>
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

    

if numel(r) > 1 && numel(minPts) > 1
    figure;
    surf(nOut)
    xticks(1:numel(minPts))
    yticks(1:numel(r))
    xticklabels(minPts)
    yticklabels(r)
    xlabel("minPts")
    ylabel("radius")
    zlabel("outlier genes")

    figure;
    surf(nClusts)
    xticklabels(minPts)
    yticklabels(r)
    xlabel("minPts")
    ylabel("radius")
    zlabel("clusters")
end

if runTopo
    profile = profile_complete(1:n,:);
    [~,OS] = generateContours2(profile(:,strcmp(sampleName,sampleClass(1))),...
        profile(:,strcmp(sampleName,sampleClass(2))),gridDensity,numContours,1);
end
    
dbOuts = numel(find(strcmp(grp,"outlier")));
topoOuts = numel(find(OS>OScutoff));
overlap = numel(intersect(find(OS>OScutoff),find(strcmp(grp,"outlier"))));

disp("-------------------------")
disp("DBSCAN unclustered genes: "+string(dbOuts))
disp("Genes w/ OS > "+string(OScutoff)+": "+string(topoOuts))
disp("-")
disp("-------------------------")
disp("Overlapping genes: "+string(overlap))

figure;
subplot(1,2,2)
p = pie([dbOuts-overlap,topoOuts-overlap,overlap]);
clr = [0.7 0.7 0.7;0 0 1;1 0 0;1 0 1];
colormap(clr(2:end,:))
labels = {'DBSCAN','AEG','both'};
legend(labels);
for i = 2 : 2 : 6
    p(i).FontSize = 15;
    p(i).FontWeight = 'bold';
end
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')


subplot(1,2,1)
tmp = intersect(find(OS>OScutoff),find(strcmp(grp,"outlier")));
grp(OS>OScutoff) = "AEG";
grp(tmp) = "both";
grpN = zeros(numel(grp),1);
for i = 1 : numel(grp)
    if strcmp(grp(i),"clustered")
        grpN(i) = 1;
    elseif strcmp(grp(i),"outlier")
        grpN(i) = 2;
    elseif strcmp(grp(i),"AEG")
        grpN(i) = 3;
    else
        grpN(i) = 4;
    end
end
lgnd = ["normal";"DBSCAN";"AEG";"both"];
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grpN,clr,'...d',[13,13,13,6])
legend(lgnd)
title(title_tmp)
xlabel(xL_tmp)
ylabel(yL_tmp)
legend('boxoff')
box off
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')


%% DBSCAN: breast cancer mouse mTOR KO treatment/control
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
[~, ind] = f.fltr(profile,0.1,3);
profile = profile(ind,:);
geneName = geneName(ind);

profile_complete = profile;
geneName_complete = geneName;




% ---------------run DBSCAN-----------------------

n = size(profile,1);


% mTOR
r = 0.7;
minPts = 80;
OScutoff = 0.1;
excludeTop = 0.05;
title_tmp = "mTOR";
xL_tmp = "Control (log_2(TPM+1))";
yL_tmp = "mTOR\Delta (log_2(TPM+1))";

runTopo = 1;
gridDensity = 100;
numContours = 20;
%--------------------------------------------------------------------------
sampleClass = ["control";"treatment"];



profileReduced = [mean(profile(:,strcmp(sampleName,sampleClass(1))),2),...
    mean(profile(:,strcmp(sampleName,sampleClass(2))),2)];

nOut = zeros(numel(r),numel(minPts));
nClusts = zeros(numel(r),numel(minPts));

for i = 1 : numel(r)
    for j = 1 : numel(minPts)
        dbOutput = dbscan(profileReduced(1:n,:),r(i),minPts(j));
        nOut(i,j) = sum(dbOutput == -1);
        nClusts(i,j) = numel(unique(dbOutput)) - 1;
        [i j]
    end
end


grp = "clust "+string(dbOutput);
grp(dbOutput == -1) = "outlier";

figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,[],[],10)
title("DBSCAN clusters")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

grp(dbOutput ~= -1) = "clustered";

figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,'cr',[],10)
title("DBSCAN classification")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

profileReduced_sorted = sort(profileReduced,1,'descend');
excludeCutoffs = profileReduced_sorted(round(n*excludeTop),:);
grp(profileReduced(:,1) > excludeCutoffs(1) &...
    profileReduced(:,2) > excludeCutoffs(2)) = "clustered";
%grp(profileReduced(:,2) > excludeCutoffs(2)) = "clustered";
figure;
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grp,'cr',[],10)
xline(excludeCutoffs(1),'m--','LineWidth',2)
yline(excludeCutoffs(2),'m--','LineWidth',2)
title("DBSCAN classification (exluding top "+excludeTop*100+"%)")
xlabel("Sample 1 mean expression (TPM)")
ylabel("Sample 2 mean expression (TPM)")
legend({"clustered","outlier"}) %#ok<CLARRSTR>
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')

    

if numel(r) > 1 && numel(minPts) > 1
    figure;
    surf(nOut)
    xticks(1:numel(minPts))
    yticks(1:numel(r))
    xticklabels(minPts)
    yticklabels(r)
    xlabel("minPts")
    ylabel("radius")
    zlabel("outlier genes")

    figure;
    surf(nClusts)
    xticklabels(minPts)
    yticklabels(r)
    xlabel("minPts")
    ylabel("radius")
    zlabel("clusters")
end

if runTopo
    profile = profile_complete(1:n,:);
    [~,OS] = generateContours2(profile(:,strcmp(sampleName,sampleClass(1))),...
        profile(:,strcmp(sampleName,sampleClass(2))),gridDensity,numContours,1);
end
    
dbOuts = numel(find(strcmp(grp,"outlier")));
topoOuts = numel(find(OS>OScutoff));
overlap = numel(intersect(find(OS>OScutoff),find(strcmp(grp,"outlier"))));

disp("-------------------------")
disp("DBSCAN unclustered genes: "+string(dbOuts))
disp("Genes w/ OS > "+string(OScutoff)+": "+string(topoOuts))
disp("-")
disp("-------------------------")
disp("Overlapping genes: "+string(overlap))

figure;
subplot(1,2,2)
p = pie([dbOuts-overlap,topoOuts-overlap,overlap]);
clr = [0.7 0.7 0.7;0 0 1;1 0 0;1 0 1];
colormap(clr(2:end,:))
labels = {'DBSCAN','AEG','both'};
legend(labels);
for i = 2 : 2 : 6
    p(i).FontSize = 15;
    p(i).FontWeight = 'bold';
end
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')


subplot(1,2,1)
tmp = intersect(find(OS>OScutoff),find(strcmp(grp,"outlier")));
grp(OS>OScutoff) = "AEG";
grp(tmp) = "both";
grpN = zeros(numel(grp),1);
for i = 1 : numel(grp)
    if strcmp(grp(i),"clustered")
        grpN(i) = 1;
    elseif strcmp(grp(i),"outlier")
        grpN(i) = 2;
    elseif strcmp(grp(i),"AEG")
        grpN(i) = 3;
    else
        grpN(i) = 4;
    end
end
lgnd = ["normal";"DBSCAN";"AEG";"both"];
gscatter(profileReduced(1:n,1),profileReduced(1:n,2),grpN,clr,'...d',[13,13,13,6])
legend(lgnd)
title(title_tmp)
xlabel(xL_tmp)
ylabel(yL_tmp)
legend('boxoff')
box off
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','bold')
