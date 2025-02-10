function run_IC_model_fmincon(putative, data)

[~, computer] = system('hostname');
if ismac
	savepath = '/Volumes/Synth-Timbre/data/manuscript/model-fits';
	% addpath('/Users/jfritzinger/Projects/WB-TIN/scripts/helper-functions',...
	% 	'-end')
elseif contains(computer, 'I1') % I1
	savepath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
else
	savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript\model-fits';
	%addpath('C:\Projects_JBF\WB-TIN\scripts\helper-functions\', '-end')
end

%% Load in data

% Load in AN
filename = sprintf('%s_AN.mat', putative);
load(fullfile(savepath, putative, filename), 'params', 'AN', 'model_params')
data_rates = analyze_data(data, CF); % Analyze data and put in correct form

%% For each AN response, fit IC model to data
timerVal2 = tic;

num_paramCF = size(AN, 1);
fit_params_all = zeros(num_paramCF,5);
for iparamCF = 1:num_paramCF
	AN_sub = AN(iparamCF,:);
	best_fval = Inf;
	for ipass = 1
		timerVal = tic;

		S1_init = 0.1 + (0.7 - 0.1) * rand(1);
		S2_init = 0.1 + (0.7 - 0.1) * rand(1);
		BMF1_init = 25 + (250 - 25) * rand(1);
		BMF2_init = 25 + (250 - 25) * rand(1);
		BMF3_init = 25 + (250 - 25) * rand(1);

		%x0 = [0.5 0.5 100 100 100]; % Slow, Shigh, D, BMFlow, BMFon, BMFhigh
		x0 = [S1_init S2_init BMF1_init BMF2_init BMF3_init];
		lb = [0   0   10  10  10 ];
		ub = [1   1   300 300 300];

		fun = @(x) minimize_IC_model_fit(...
			data_rates, AN_sub, params, model_params, x);
		options = optimoptions('fmincon','Display','off',...
			'Algorithm','sqp');
		[fit_params, fval] = fmincon(...
			fun, x0, [], [], [], [], lb, ub, [], options);
		

		if fval < best_fval
			best_x = fit_params;
			best_fval = fval;
		end

		fprintf('Fitting took %0.02f minutes\n', toc(timerVal)/60);
	end
	fit_params_all(iparamCF, :) = best_x;
end

% Save IC model fit results
filename = sprintf('%s_IC.mat', putative);
save(fullfile(savepath, putative, filename), 'fit_params_all')
fprintf('Fitting all models took %0.02f minutes\n', toc(timerVal2)/60);


