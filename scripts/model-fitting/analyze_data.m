function data_rates = analyze_data(data, CF)

%% For first fmincon try 'AN' and 'AN_morereps'
% params_WB = data(7,2); % Gets WB-TIN 23 dB 
% params_NB = data(10, 2); % Gets NB-TIN 23 dB
% params_MTF = data{3,2}; % Gets binaural MTFN
% data_MTF = analyzeMTF(params_MTF);
% data_WB = analyzeWBTIN(params_WB, []);
% data_WB = data_WB{1};
% data_NB = analyzeNBTIN(params_NB, CF);
% data_NB = data_NB{1};
% 
% rate_MTF = data_MTF.rate_sm; % smoothed
% rate_TIN = data_NB.rate_onCF'; % not smoothed 
% rate_WB = data_WB.rate(:,2); % not smoothed
% 
% data_rates = [rate_MTF; rate_TIN; rate_WB];

%% For fmincon to fit to synthetic timbre data 

params_WB = data(7,2); % Gets WB-TIN 23 dB 
params_MTF = data{3,2}; % Gets binaural MTFN
data_MTF = analyzeMTF(params_MTF);
data_WB = analyzeWBTIN(params_WB, []);
data_WB = data_WB{1};

rate_MTF = data_MTF.rate_sm; % smoothed
rate_WB = data_WB.rate(:,2); % not smoothed

data_rates = [rate_MTF; rate_WB];

end