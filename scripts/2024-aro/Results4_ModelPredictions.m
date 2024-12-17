%% Results4_ModelPredictions
% J. Fritzinger, 1/16/23
clear

%% Phyiso session matched to model

% BS Example 1
% session = 'R024S470_TT2_N1';
% CF = 1200;
% iMTF = 2;

figure('Position', [313,693,1263,513])
tiledlayout(2,3 , 'TileSpacing','compact', 'Padding','compact')
fontsize = 16;

for ineuron = 1:2

	switch ineuron
		case 1
			% BE Example 1
			% R027S247, TT3N1, Excellent, CF = 1741Hz, BE
			session = 'R027S247_TT3_N1';
			CF = 1500;
			iMTF = 1;
		case 2
			% R029S105, TT1N1, Good, CF = 6063Hz, BE
			session = 'R029S105_TT1_N1';
			CF = 6063;
			iMTF = 1;
	end

	%% Set up paths

	% Load in data
	base = getPaths();
	fpath = 'data/aro-2024';
	filename = fullfile(base, fpath, [session '.mat']);
	load(filename, 'params', 'data', 'cluster', 'stims')

	% Synth Timbre
	has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,...
		'Spectral_Centroid')&&p.binmode==2&&p.Delta_F==200&&mod(p.fpeak_mid, 200)==0,params);
	data = data(has_sc); % Gets binaural WB-TIN stimuli
	params = params(has_sc);
	num_ds = size(params, 1);

	%% Set up synth timbre stimulus
	% for ind = 1:size(params, 1)
	% 	params{ind}.Fs = 100000;
	% 	params{ind}.physio = 1;
	% 	params{ind}.mnrep = 30;
	% 	params{ind}.dur = 0.3;
	% 	[params{ind}] = generate_synthetictimbre(params{ind});
	% 	params{ind}.num_stim = size(params{ind}.stim, 1);
	% end
	% 
	% % Set up MTF
	% params{ind+1}.Fs = 100000;
	% params{ind+1}.mnrep = 5;
	% params{ind+1}.physio = 0;
	% params{ind+1}.type = 'typMTFN';
	% params{ind+1}.ramp_dur = 0.05;
	% params{ind+1}.noise_state = 0;
	% params{ind+1}.noise_band = [100, 10000];
	% params{ind+1}.dur = 1;
	% params{ind+1}.reptim = 1.5;
	% params{ind+1}.fms = [2, 600, 3];
	% params{ind+1}.mdepths = [0,0,1];
	% params{ind+1}.binmode = 2;
	% params{ind+1}.No = 30;
	% params{ind+1}.spl = 30;
	% params{ind+1}.raised_sine = 1;
	% params{ind+1}.onsetWin = 25;
	% params{ind+1} = generate_MTF(params{ind+1});
	% params{ind+1} = generate_MTF(params{ind+1});
	% params{ind+1}.num_stim = size(params{ind+1}.stim, 1);
	
	% %% Run model
	% 
	% % Set up cell arrays for model results
	% num_stim = size(params, 1);
	% energy = cell(num_stim, 1);
	% lateral_model = cell(num_stim, 1);
	% for ist = 1:num_stim
	% 
	% 	% Energy model
	% 	timerVal = tic;
	% 
	% 	CFs = CF;
	% 	Fs = params{ist}.Fs;
	% 	stimulus = [params{ist}.stim zeros(size(params{ist}.stim,1),0.1*Fs)];
	% 	gamma_param.srate = Fs;
	% 	tvals = (1:length(stimulus))/Fs;
	% 	gamma_IF_reg = zeros(length(CFs),length(tvals));
	% 	impaired = 0; % 0 = not impaired; 1 = 'impaired'
	% 
	% 	pin_gamma = zeros(size(stimulus, 1), Fs*params{ist}.dur+0.1*Fs);
	% 	for istim = 1:size(stimulus, 1)
		% 	gamma_param.fc = CFs;
		% 	pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
	% 	end
	% 	pin_gamma = pin_gamma(:,1:params{ist}.dur*Fs);
	% 	energy{ist} = sqrt(mean(pin_gamma.^2,2));
	% 
	% 	elapsedTime = toc(timerVal)/60;
	% 	disp(['Energy model took ' num2str(elapsedTime) ' minutes'])
	% 
	% 	% Set up lateral inhibition model and run
	% 	timerVal = tic;
	% 	clear model_params
	% 
	% 	% Lateral Model
	% 	model_params.type = 'Lateral Model';
	% 	model_params.lateral_CF = [CF/2 CF CF*2];
	% 	model_params.CFs = model_params.lateral_CF;
	% 	model_params.CF_range = model_params.CFs(2);
	% 	if iMTF == 1 % BE
		% 	CS_param = [0.45, 0.45, 0.003];
		% 	model_params.config_type = 'BS inhibited by off-CF BS';
	% 	else % BS
		% 	CS_param = [0.08, 0.08, 0.003];
		% 	model_params.config_type = 'BS inhibited by off-CF CN';
	% 	end
	% 
	% 	% Model parameters
	% 	model_params.range = 2; % 1 = population model, 2 = single cell model
	% 	model_params.species = 1; % 1 = cat, 2 = human
	% 	model_params.num_CFs = 1;
	% 	model_params.nAN_fibers_per_CF = 5;
	% 	model_params.cohc = 1; % (0-1 where 1 is normal)
	% 	model_params.cihc = 1; % (0-1 where 1 is normal)
	% 	model_params.nrep = 1; % how many times to run the AN model
	% 	model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
	% 	model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
	% 	model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
	% 	model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
	% 	model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
	% 	model_params.BMF = 100;
	% 
	% 	% Run Lateral Inhibition Model
	% 	AN = modelLateralAN(params{ist}, model_params);
	% 	lateral_model{ist} = modelLateralSFIE(params{ist}, ...
		% 	model_params, AN.an_sout, AN.an_sout_lo, AN.an_sout_hi,'CS_params', CS_param);
	% 
	% 	elapsedTime = toc(timerVal)/60;
	% 	disp(['Lateral inhibition model took ' num2str(elapsedTime) ' minutes'])
	% end
	% 
	% %% Save model results
	% savename = fullfile(path, [session '_LatModel.mat']);
	% save(savename, 'params', 'energy', 'model_params', 'lateral_model', 'data', '-v7.3')

	%% Load in data

	savename = fullfile(base, fpath, [session '_LatModel.mat']);
	load(savename)

	%% Plot
	if num_ds == 4
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};
	else
		data_colors = {'#82BB95', '#03882F', '#034E1C'};
		model_colors = {'#E29A63', '#984102', '#4E2201'};
	end
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];

	ylimits = [0 80];

	% Plot model MTF
	[rate_MTF, rate_std] = plotMTF(params{4}, lateral_model{4}.avIC, 0);

	nexttile
	hold on
	line([1 params{4}.all_fms(end)], [1 1]*rate_MTF(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
	errorbar(params{4}.all_fms,rate_MTF,rate_std,'Color','k', 'Marker','.', 'MarkerSize',10, 'LineWidth', 2);

	% Label the plots
	box on
	grid on
	xtick = [1 2 5 10 20 50 100 200 500];
	xlim(xtick([1 end]))
	xlabel('Modulation Freq (Hz)')
	ylabel('Spike Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	hold off
	set(gca,'fontsize',fontsize)
	legend('Unmodulated', 'fontsize',fontsize, 'location', 'northwest', 'EdgeColor','none')
	if ineuron ==1
	title('Lateral Model MTF', 'FontSize', 18);
	end

	% Plot lateral inhibition
	clear leg
	nexttile
	hold on
	for ind = 1:3
		[rate, ~] = plotSyntheticTimbre(params{ind}, lateral_model{ind}.avIC, 0);
		plot(params{ind}.fpeaks,rate, ...
			'linewidth', 1.5, 'Color', model_colors{ind});
		R_int = corrcoef(data{ind}.rate,rate);
		R(ind) = R_int(1, 2).^2;
		leg{ind} = [num2str(params{ind}.spl) ' dB SPL, R^2 = ' num2str(round(R(ind), 3))];
	end
	leg{ind+1} = 'CF';
	xlabel('Spectral Peak Frequency (Hz)')
	ylim(ylimits)
	set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	%yticklabels([])
	ylabel('Avg. Rate (sp/s)')
	xlim(plot_range);
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
	%legend(leg,'location', 'northeast', 'FontSize',8)

	set(gca,'FontSize',fontsize)
	box on
	grid on
	if ineuron == 1
	title('Lateral Inhibition Model', 'FontSize',18)
	end

		% R^2 annotation
	for ind = 1:num_ds
		if ineuron == 1
			start = 0.85;
		else
			start = 0.34;
		end
		annotation('textbox',[0.57 start-0.033*(ind-1) 0.2 0.0917],...
			'Color',model_colors{ind},...
			'String',['R^2 = ' num2str(round(R(ind), 2))],...
			'FontSize',16,...
			'FitBoxToText','off',...
			'EdgeColor','none');
		annotation('textbox',[0.37 start-0.033*(ind-1) 0.2 0.0869],...
			'Color',model_colors{ind},...
			'String',[num2str(params{ind}.spl) ' dB SPL'],...
			'FontSize',16,...
			'EdgeColor','none');
	end

	% Plot data
	clear leg
	nexttile
	hold on
	for ind = 1:3
		errorbar(data{ind}.fpeaks,data{ind}.rate,data{ind}.rate_std/(sqrt(params{ind}.nrep)),...
			'LineWidth', 1, 'Color', data_colors{ind});
		leg(ind,:) = [num2str(params{ind}.spl) ' dB SPL'];
	end
	%plot(data{ind}.fpeaks{1}([1 end]),[1 1]*data{ind}.spont,'-','Color',[0.5 0.5 0.5], 'LineWidth', 2)
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
	%legend(char(leg, 'CF'),...
	%	'location', 'northeast', 'FontSize',12)
	plot_range = [data{ind}.fpeaks(1) data{ind}.fpeaks(end)];
	xlabel('Spectral Peak Frequency (Hz)')
	xlim(plot_range);
	set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	ylabel('Avg. Rate (sp/s)')
	ylim(ylimits)
	set(gca,'FontSize',fontsize)
	box on
	grid on
	if ineuron == 1
	title('Single-Unit', 'FontSize',18);
	end
	
	for ind = 1:num_ds
		annotation('textbox',[0.69 start-0.033*(ind-1) 0.2 0.0869],...
			'Color',data_colors{ind},...
			'String',[num2str(params{ind}.spl) ' dB SPL'],...
			'FontSize',16,...
			'EdgeColor','none');
	end

	% % Energy
	% clear leg
	% nexttile
	% hold on
	% for ind = 1:3
	% 	[rate, ~] = plotSyntheticTimbre(params{ind}, energy{ind}, 0);
	% 	plot(params{ind}.fpeaks,rate, ...
		% 	'linewidth', 1.5, 'Color', model_colors{ind});
	% 	R_int = corrcoef(data{ind}.rate,rate);
	% 	R(ind) = R_int(1, 2).^2;
	% 	leg(ind,:) = [num2str(params{ind}.spl+3) ' dB SPL, R^2 = ' num2str(round(R(ind), 3))];
	% end
	% xlabel('Spectral Peak Frequency (Hz)')
	% %ylim([0 0.003])
	% set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	% xlim(plot_range);
	% xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
	% legend(char(leg,'CF'),'location', 'northeast', 'FontSize',12)
	% set(gca,'FontSize',14)
	% box on
	% grid on
	% title('Energy', 'FontSize',18)
	% ylabel('Avg. Rate (sp/s)')


end







