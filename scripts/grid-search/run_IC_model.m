function run_IC_model(putative_timbre)

% Set up grid search parameters
paramS_range = 0:0.05:0.6;
paramS2_range = [];
paramCF_range = 0.25:0.25:1.25;
num_paramS = length(paramS_range);
num_paramS2 = length(paramS2_range);
num_paramCF = length(paramCF_range);
paramD = 0;

% Run model
timerVal = tic;
lateral_model = cell(num_paramS, 3);

for iparamCF = 1:num_paramCF
	timerVal5 = tic;

	% Load AN model
	filename = sprintf('%s_AN_%d.mat', putative_timbre, iparamCF);
	load(fullfile(savepath, putative_timbre, filename), 'params', 'AN', 'CF', 'model_params')

	for iparamS = 1:num_paramS  % Run IC model
		lm_params = [paramS_range(iparamS) paramS_range(iparamS) paramD];
		parfor ist = 1:num_stim
			if ~isempty(params{ist})
				lateral_model{iparamS, ist} = modelLateralSFIE(params{ist}, ...
					model_params, AN{ist}.an_sout, AN{ist}.an_sout_lo, ...
					AN{ist}.an_sout_hi,'CS_params', lm_params);
				lateral_model{iparamS, ist}.ic = [];
			end
		end
	end

	filename = sprintf('%s_IC_%d.mat', putative_timbre, iparamCF);
	save(fullfile(savepath, putative_timbre, filename), ...
		'params', 'model_params', 'lateral_model', 'paramS_range',  ...
		'paramCF_range', '-v7.3')
	disp(['IC Model ' num2str(iparamCF) ' took ' num2str(toc(timerVal5)) ' seconds'])

end
disp(['All IC models took ' num2str(toc(timerVal)/60) ' minutes'])

end