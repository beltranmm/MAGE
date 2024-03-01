
%% setup

% load MAGE functions
f = MAGE_functions;

profile = f.simData(1000,4,7,1,0.1);

figure;
scatter(mean(profile(:,1:4),2),mean(profile(:,5:8),2))

% filter code here

%% run MAGE

[~,OS] = MAGE(profile(:,1:4),profile(:,5:8),100,5,1);