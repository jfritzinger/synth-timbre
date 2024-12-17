%% exploration_GLM_models.m
clear

%% Load in example dataset 

[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Synth-Timbre/data/manuscript';
putative = 'R24_TT2_P13_N05';

load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')

%% Analyze/set up predictors 

% Predictors 
% Column = one predictor variable
% Row = one observation
% Load in model data 
X(:,1) = energy{2,1}.rate;
X(:,2) = SFIE{2,1}.rate;



% Response variable 
% Column vector with same number of rows as X 
load(fullfile(datapath, 'neural_data', [putative '.mat']))
params_ST = data(7, 2);
data_ST = analyzeST(params_ST);
data_ST = data_ST{1};

y = data_ST.rate;


%% GLM 

mdl = glmfit(X, y);
disp(mdl)
yfit = glmval(mdl,X, 'identity');

%% 

figure
plot(X, y, 'o')
hold on
plot(X, yfit)