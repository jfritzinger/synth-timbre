%% exploration_GLM_neural.m
clear

%% Load in example dataset 

[base, datapath, savepath, ppi] = getPaths();
putative = 'R29_TT3_P5_N10';

%% Analyze/set up predictors 

% Predictors 
% Column = one predictor variable
% Row = one observation

% STRF prediction
STRF_dataset = 'R29_TT3_P5_N10_STRF_63_Bin.mat';
load(fullfile(datapath, 'STRF_Models', STRF_dataset));


X = STRFmodel.avModel;

% Response variable 
% Column vector with same number of rows as X 
load(fullfile(datapath, 'neural_data', [putative '.mat']))
params_ST = data(7, 2);
data_ST = analyzeST(params_ST);
data_ST = data_ST{1};

y = data_ST.rate;

%% GLM 

[mdl,dev,stats] = glmfit(X, y, 'normal');
disp(mdl)


%% Analyze output 

yfit = glmval(mdl,X, 'identity');

figure
plot(X, y, 'o')
hold on
plot(X, yfit)

figure
plot(data_ST.fpeaks, data_ST.rate)
hold on
plot(data_ST.fpeaks, STRFmodel.avModel.*STRFmodel.ratio)