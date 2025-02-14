function get_best_fit_model(putative, CF, putative_timbre)

[~, computer] = system('hostname');
if ismac
	modelpath = '/Volumes/WB-TIN/data/model-fits';
	addpath('/Users/jfritzinger/Projects/WB-TIN/scripts/helper-functions',...
		'-end')
elseif contains(computer, 'I1') % I1
    modelpath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
    addpath('\\NSC-LCARNEY-H2\Projects_JBF\WB-TIN\scripts\helper-functions\', '-end')
else
	modelpath = 'C:\DataFiles_JBF\WB-TIN\data\model-fits';
	addpath('C:\Projects_JBF\WB-TIN\scripts\helper-functions\', '-end')
end
[~, datapath, ~, ~] = getPathsWBTIN();

%% Load in data 

load(fullfile(datapath, [putative '.mat']), 'data');
data_rates = analyze_data(data, CF); % Analyze data and put in correct form 


%% Load in AN response

% Load in AN
filename = sprintf('%s_AN.mat', putative_timbre);
load(fullfile(modelpath, putative_timbre, filename), 'params', 'AN', 'model_params')

% Load in IC parameter values 
filename = sprintf('%s_IC.mat', putative_timbre);
load(fullfile(modelpath, putative_timbre, filename), 'fit_params_all')

%% Evaluate model fits  

data_MTF = data_rates(1:26);
data_TIN = data_rates(27:29);
data_WB = data_rates(30:end);
num_paramCF = size(AN, 1);
for iparamCF = 1:num_paramCF

	% Run model with fit parameters
	AN_sub = AN(iparamCF,:);
	fit_params = fit_params_all(iparamCF,:);
	CS_params = [fit_params(1:2) 0.001];
	BMFs = fit_params(3:5);
	nstim = size(params, 2);
	model_outputs = cell(nstim, 1);
	for istim = 1:nstim
		param = params{istim};
		an_sout = squeeze(AN_sub{istim}.an_sout);
		an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
		an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
		model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
			an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
			'BMFs', BMFs);
	end

	% Analyze model output
	for istim = 1:nstim
		param = params{istim};
		model_output = model_outputs{istim};
		switch param.type
			case 'typMTFN'
				[~, model_MTF, ~, ~] = plotMTF(param, model_output.avIC, 0);
			case 'TIN'
				[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
			case 'SPEC_slide' % WB-TIN only for now
				[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
		end
	end

	% Calculate MSE 
	mse = minimize_IC_model_fit(data_rates, AN_sub, params, model_params,...
		fit_params);
	mse_all(iparamCF) = mse;

end

% Get best parameters
[mse_best, min_ind] = min(mse_all);
best_CF_range = AN{min_ind, 1}.CF_span;
best_params = 1;

%% Run model with best parameters
 
paramCF = min_ind;
AN_sub = AN(paramCF,:);
fit_params = fit_params_all(paramCF,:);
CS_params = [fit_params(1:2) 0.001];
BMFs = fit_params(3:5);
nstim = size(params, 2);
IC_best = cell(nstim, 1);
for istim = 1:nstim
	param = params{istim};
	an_sout = squeeze(AN_sub{istim}.an_sout);
	an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
	an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
	IC_best{istim} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
end
AN_best = AN(paramCF, :);

%% Save model and best parameters

filename = sprintf('%s_BestModel.mat', putative);
save(fullfile(modelpath, putative_timbre, filename),...
	'AN_best', "IC_best", "model_params", "params", "fit_params")

