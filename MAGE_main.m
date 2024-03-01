
%% setup

% load MAGE functions
f = MAGE_functions;

profile = f.simData(10000,4,7,1,0.1);

figure;
scatter(mean(profile(:,1:4),2),mean(profile(:,5:8),2))

filtered_profile = f.fltr(profile,5,4);
figure;
scatter(mean(filtered_profile(:,1:4),2),mean(filtered_profile(:,5:8),2))