% plot_model_intermediate
clear

% param.type = 'MTFN';
% param.Fs = 100000;
% param.mnrep = 10;
% param.ramp_dur = 0.05;
% param.noise_state = 0;
% param.noise_band = [100, 10000];
% param.dur = 1;
% param.reptim = 1.5;
% param.fms = [2, 600, 2];
% param.mdepths = [0,0,1];
% param.binmode = 2;
% param.No = 30;
% param.spl = 30;
% param.raised_sine = 1;
% param.onsetWin = 25;
% param = generate_MTF(param);
% param.num_stim = size(param.stim, 1);

%% Run model

IC = cell(2,1);
% range = 0.75;
% CS_params = [0.5 0.5 0]; %0.0039];
% BMFs = [100 100 100];
for iexample = 1:2

	switch iexample
		case 1
			range = 0.62;
			CS_params = [1 0.28 0];
			BMFs = [87 164 111];
			CF = 1000;
		case 2
			range = 0.5;
			CS_params = [0.19 0.00 0.001];
			BMFs = [300	38	10];
			CF = 3482;
			%CF = 1000;
	end

	% Generate SAM tone or noise
	param.type = 'MTFT';
	param.Fs = 100000;
	param.ramp_dur = 0.05;             % All time in seconds
	param.carrier_freq = CF;
	param.spl = 70;
	param.binmode = 1;
	param.w_bgn = 0;
	param.No = 20;
	param.nrep = 3; %10
	param.dur = 1;                 % Stimulus duration
	param.reptim = 1.5;            % Period between stimulus onsets
	param.fms = [2 600 2]; % [4 1025 3];        % lo, hi, and steps/oct.
	param.mdepths = [0 0 1];       % lo, hi, and steps in dB.
	param = generate_MTFT(param);
	param.num_stim = size(param.stim, 1);

	% Model parameters
	model_params.range = 2; % 1 = population model, 2 = single cell model
	model_params.species = 1; % 1 = cat, 2 = human
	model_params.num_CFs = 1;
	model_params.nAN_fibers_per_CF = 5;
	model_params.cohc = 1; % (0-1 where 1 is normal)
	model_params.cihc = 1; % (0-1 where 1 is normal)
	model_params.nrep = 1; % how many times to run the AN model
	model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
	model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
	model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
	model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
	model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
	model_params.BMF = 100;
	model_params.type = 'Lateral Model';
	model_params.config_type = 'BS inhibited by off-CF BS';

	% Run AN Model
	timerVal = tic;
	model_params.lateral_CF = [CF*2^(-1*range), CF, CF*2^range];
	model_params.CFs = model_params.lateral_CF;
	model_params.CF_range = model_params.CFs(2);

	% Used for manuscript
	AN = modelLateralAN(param, model_params);
	disp(['AN model took ' num2str(toc(timerVal)) ' seconds'])

	% Run IC model
	an_sout = squeeze(AN.an_sout);
	an_sout_lo = squeeze(AN.an_sout_lo);
	an_sout_hi = squeeze(AN.an_sout_hi);
	IC{iexample} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
	IC{iexample}.CFs = [CF*2^(-1*range), CF, CF*2^range];
	IC{iexample}.CS_params = CS_params;
	IC{iexample}.BMFs = BMFs;
end

% Save model responses
% [base, datapath, ~, ~] = getPathsWBTIN();
% filename = 'Model_Component_Examples.mat';
% save(fullfile(base, filename), 'param', 'IC', ...
% 	'model_params', '-v7.3')
% 

%% Analysis and Plotting

fig = figure('Position',[45,400,947,524]);
%tiledlayout(2, 4)
fontsize = 14;
titlesize = 16;
labelsize = 22;

index = [1, 2, 3, 4; 5, 6, 7, 8];
for iexample = 1:2

	[~, model_MTF(:,1), model_std(:,1), ~, smoothed(:,1)] = plotMTF(param, IC{iexample}.avBS_lo, 0);
	[~, model_MTF(:,2), model_std(:,2), ~, smoothed(:,2)] = plotMTF(param, IC{iexample}.avBS_on, 0);
	[~, model_MTF(:,3), model_std(:,3), ~, smoothed(:,3)] = plotMTF(param, IC{iexample}.avBS_hi, 0);
	[~, model_MTF(:,4), model_std(:,4), ~, smoothed(:,4)] = plotMTF(param, IC{iexample}.avIC, 0);
	CFs = IC{iexample}.CFs;
	CS_params = IC{iexample}.CS_params;

	% Plot MTF
	names = {'Low BS', 'On BS', 'High BS', 'Broad Inhibition', 'On BE'};
	for ii = 1:3
		h(iexample, ii) = subplot(3, 4, index(iexample, ii));
		hold on
		yline(model_MTF(1,ii),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
		param.all_fms(1) = 1.2;
		errorbar(param.all_fms, model_MTF(:,ii), model_std(:,ii), 'linewidth', 1.5)
		hold off
		xtick = [1 2 5 20 50 100 200 500];
		xlim([1 530])
		xlabel('Modulation Freq (Hz)')
		if ii == 1
			ylabel('Avg. Rate (sp/s)')
		else
			yticklabels([])
		end
		set(gca,'XTick',xtick,'XScale', 'log')
		if iexample == 2 && ii == 1
			legend('Unmodulated', 'Location','best')
		end
		grid on
		axis_set = axis;
		if iexample == 1
			axis_set(3) = 0;
			axis_set(4) = 50;
		else
			axis_set(3) = 0;
			axis_set(4) = 60;
		end

		if ii == 2
			if iexample == 1
				CS_params = [1 0.28 0];
			else
				CS_params = [0.19 0.00 0.001];
			end
			modifiers = CS_params(1:2);

			hold on
			plot(param.all_fms, model_MTF(:,1).*modifiers(1), '--')
			plot(param.all_fms, model_MTF(:,3)*modifiers(2), ':')
			plot(param.all_fms, model_MTF(:,4), 'k')
		end

		axis(axis_set);
		set(gca, 'FontSize', fontsize)
		if ii == 4 || ii == 5
			title(names{ii}, 'FontSize',titlesize)
		else
			title([names{ii} ', CF = ' num2str(round(CFs(ii))) 'Hz'], 'FontSize',titlesize)
		end
	end

	h(iexample, 4) = subplot(3, 4, index(iexample, 4));
	rate_MTF = model_MTF(:,2) - CS_params(1)*model_MTF(:,1) - CS_params(2)*model_MTF(:,3);
	hold on
	yline(model_MTF(1,4),'Color','k', 'linewidth', 1.5);
	errorbar(param.all_fms, model_MTF(:,4), model_std(:,4),'k', 'linewidth', 1.5)
	%yline(rate_MTF(1),'Color','r', 'linewidth', 1.5);
	%plot(param.all_fms,rate_MTF,'r', 'linewidth', 1.5);
	hold off
	xtick = [1 2 5 20 50 100 200 500];
	xlim([1 530])
	xlabel('Modulation Freq (Hz)')
	ylabel('Avg. Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	%legend('Unmodulated', 'Location','best')
	grid on
	set(gca, 'FontSize', fontsize)
	% axis_set = axis;
	% axis_set(3) = 0;
	% axis_set(4) = 30;
	% axis(axis_set);
	if iexample == 2
		ylim([0 60])
	end
	if iexample == 1
		legend('Broad Inhibition MTF', '', 'Rate Subtracted MTF', '')
	end
	title('Comparison', 'FontSize',titlesize)

end

%% Load in images
% 
% img = imread('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/WB-TIN/data/2024-manuscript/Fig14_Data/BroadInhModel_Split.png');
% 
% p(1) = subplot(3, 4, 9);
% imshow(img(1:566,2:687,:))
% 
% p(2) = subplot(3, 4, 10);
% imshow(img(567:1134,2:687,:))
% 
% p(3) = subplot(3, 4, 11);
% imshow(img(1135:end,2:687,:))

%% Arrange 

left = [0.08 0.28 0.48 0.75];
bottom = [0.62 0.11];
width = 0.18;
height = 0.3;

for iexample = 1:2
	for ii = 1:4
		set(h(iexample,ii), 'Position', [left(ii) bottom(iexample) width height])
	end
end

set(p(1), 'Position',[0.085 0.62 0.13 0.19])
set(p(2), 'Position',[0.285 0.62 0.13 0.19])
set(p(3), 'Position',[0.485 0.62 0.125 0.19])

annotation('textbox',[left(1)-0.05 0.95 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(4)-0.05 0.95 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1)-0.05 0.44 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(4)-0.05 0.44 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Example Neuron
%
% [~, computer] = system('hostname');
% if ismac
% 	modelpath = '/Volumes/WB-TIN/data/model-fits';
% 	addpath('/Users/jfritzinger/Projects/WB-TIN/scripts/helper-functions',...
% 		'-end')
% elseif contains(computer, 'I1') % I1
% 	modelpath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\model-fits';
% else
% 	modelpath = 'C:\DataFiles_JBF\WB-TIN\data\model-fits';
% 	addpath('C:\Projects_JBF\WB-TIN\scripts\helper-functions\', '-end')
% end
% [base, datapath, ~, ~] = getPathsWBTIN();
%
% %% Load in data
%
% name = 'BE';
% putative = 'R29_TT4_P2_N16';
% CF = 2639;
%
% load(fullfile(datapath, [putative '.mat']), 'data');
% data_rates = analyze_data(data, CF); % Analyze data and put in correct form
%
% % Load in AN
% filename = sprintf('%s_AN_old.mat', putative);
% load(fullfile(modelpath, putative, filename), 'param', 'AN', 'model_param')
%
% % Load in IC parameter values
% filename = sprintf('%s_IC_old.mat', putative);
% load(fullfile(modelpath, putative, filename), 'fit_param_all')
%
% %% Get best parameters
%
% iparamCF = 7;
%
% % Run model with fit parameters
% AN_sub = AN(iparamCF,:);
% fit_param = fit_param_all(iparamCF,:);
% CS_param = [fit_param(1:2) 0.001];
% %BMFs = fit_param(3:5);
% BMFs = fit_param(4:6);
% nstim = size(param, 2);
% model_outputs = cell(nstim, 1);
%
% timerVal2 = tic;
% for istim = 1 %1:nstim
% 	param = param{istim};
% 	an_sout = squeeze(AN_sub{istim}.an_sout);
% 	an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
% 	an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
% 	model_outputs{istim} = modelLateralSFIE_BMF(param, model_param, ...
% 		an_sout, an_sout_lo, an_sout_hi, 'CS_param', CS_param,...
% 		'BMFs', BMFs);
% end
% disp(['IC model took ' num2str(toc(timerVal2)) ' seconds'])
%
% %% Analyze model outputs, intermediate and final steps
%
% for istim = 1 % 1:nstim
% 	param = param{istim};
% 	model_output = model_outputs{istim};
% 	switch param.type
% 		case 'typMTFN'
% 			[~, model_MTF(:,1), ~, ~] = plotMTF(param, model_output.avBS_lo, 0);
% 			[~, model_MTF(:,2), ~, ~] = plotMTF(param, model_output.avBS_on, 0);
% 			[~, model_MTF(:,3), ~, ~] = plotMTF(param, model_output.avBS_hi, 0);
% 			[~, model_MTF(:,4), ~, ~] = plotMTF(param, model_output.avIC, 0);
% 		case 'TIN'
% 			[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
% 		case 'SPEC_slide' % WB-TIN only for now
% 			[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
% 			model_WB = model_WB(:,2);
% 	end
% end
%
% %% Plot
%
% fig = figure();
% tiledlayout(1, 4, 'Padding','compact')
%
% % Plot MTF
% for ii = 1:4
% 	nexttile
% 	param = param{1};
% 	hold on
% 	%yline(data_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
% 	%plot(param.all_fms,data_MTF);
% 	yline(model_MTF(1,ii),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
% 	plot(param.all_fms,model_MTF(:,ii));
% 	%plot(data.fms,rate_sm,'-b', 'LineWidth', 1)
% 	hold off
% 	xtick = [1 2 5 20 50 100 200 500];
% 	xlim(xtick([1 end]))
% 	xlabel('Modulation Freq (Hz)')
% 	ylabel('Avg. Rate (sp/s)')
% 	set(gca,'XTick',xtick,'XScale', 'log')
% 	legend('Unmodulated', 'Location','best')
% 	grid on
% 	axis_set = axis;
% 	axis_set(3) = 0;
% 	axis(axis_set);
% end

