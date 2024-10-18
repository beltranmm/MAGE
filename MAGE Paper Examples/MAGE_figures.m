%% setup

% select data directory
dataDir = strcat(pwd,'\workspaces\');

% load contour functions
f = DEG_contour_functions;


profile_complete = zeros();

%% figure 3 [DBSCAN]


% load figure data

% generate figure

%% figure 4 [MAGE g-T3]


% load figure data

% generate figure

%% figure 5 [Comparison g-T3]


% load figure data
fig = openfig(strcat(pwd,'\figures\gt3 comparison.fig'));
axObjs = fig.Children
dataObjs = axObjs.Children

x = dataObjs(3).XData
y = dataObjs(3).YData

% generate figure

%% figure 6 [MAGE mTOR]


% load figure data

% generate figure

%% figure 7 [Comparison mTOR]


% load figure data

% generate figure

