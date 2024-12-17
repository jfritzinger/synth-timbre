% Results1_BE_Example
% J. Fritzinger, 1/16/23
clear

%% Plot datasets

for ineuron = 1:2

	switch ineuron
		case 1

			% R027S247, TT3N1, Excellent, CF = 1741Hz, BE
			session = 'R027S247_TT3_N1';
			CF = 1500;
			ylimits = [0 60];

			% R027S037, TT3N1, Good, CF = 1741Hz, BE
			% session = 'R027S037_TT3_N1';
			% CF = 1741;
			% ylimits = [0 60];

			% R027S292, TT3N1, Excellent, CF = 1741Hz, BE
			% session = 'R027S292_TT3_N1';
			% CF = 1741;

		case 2
			
			% R025S467, TT2N1, Excellent, CF = 2997Hz
			% session = 'R025S467_TT2_N1';
			% CF = 2997;

			% % R025S474, TT1N3, Good, CF = 4595Hz, BE 
			% session = 'R025S474_TT1_N3';
			% CF = 4595;
			% 
			% R029S045, TT3N1, Good, CF = 8000Hz, BE
			% session = 'R029S045_TT3_N1';
			% CF = 8000;
			
			% R029S099, TT3N1, Good, CF = 6063Hz, BE
			% session = 'R029S099_TT3_N1';
			% CF = 6063;

			% R029S105, TT1N1, Good, CF = 6063Hz, BE
			session = 'R029S105_TT1_N1';
			CF = 6063;
			ylimits = [0 55];

	end

	%figure('Position', [282,533,1474,345]);
	figure('Position',[288,108,1250,360])
	t = tiledlayout(4, 2, 'TileSpacing','compact', 'TileIndexing','columnmajor');
	spont_color = [0.4 0.4 0.4];
	CF_color = [0.3 0.3 0.3];
	x_label = [200 500 1000 2000 5000 8000];

	%% Load in examples

	base = getPaths();
	fpath = 'data/aro-2024';
	filename = fullfile(base, fpath, [session '.mat']);
	load(filename, 'params', 'data', 'cluster', 'stims')

	% RM
	data_RM = data{4};

	% MTF
	params_MTF = params{5};
	data_MTF = data{5};

	% STRF
	%has_strf = cellfun(@(p) strcmp(p.type,'STRF'),params);
	%data_STRF = data{has_strf};

	% Synth Timbre
	has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,...
		'Spectral_Centroid')&&p.binmode==2&&p.Delta_F==200&&mod(p.fpeak_mid, 200)==0,params);
	data_ST = data(has_sc); % Gets binaural WB-TIN stimuli
	params_ST = params(has_sc);
	plot_range = [params_ST{1}.fpeaks(1) params_ST{1}.fpeaks(end)];


	%% Response map
	color_RM = {'#A49BD0', '#5E50A9', '#20116B'};


	% Plot
	%nexttile(t)
	for ind = 1:3
		h(ind) = subplot(1, 6, ind);
		hold on
		plot(data_RM.freqs,data_RM.rates(:,ind+2),'color', color_RM{ind},'LineWidth',2) % 33 dB
		plot(data_RM.freqs([1 end]),[1 1]*data_RM.spont,'-','LineWidth',2, 'Color',spont_color)
		xline(CF, '--', 'Color', CF_color,'LineWidth',2);
		box on
		grid on
		hold off
		ylim([0 max(data_RM.rates, [], 'all')+20])
		set(gca,'XTick',[])
		xlim([250 plot_range(2)])
		xticks(x_label)
		set(gca, 'XScale', 'log');
		set(gcf, 'color', 'w')
		set(gca,'fontsize',14)

		if ind == 1
			xlabel('Tone Frequency (Hz)')
		elseif ind == 2
			xticklabels([])
			ylabel('Avg. Rate (sp/s)')
		elseif ind == 3
			if ineuron == 1
				title('Response Map', 'FontSize', 18)
			end
			xticklabels([])
		end

	end

	%% MTF
	MTF_range = [0 max(data_MTF.rate)+5];

	% Plot
	%nexttile(t)
	h(4) = subplot(1, 6, 4);
	hold on
	line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',2);
	errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',2, 'Color','k');
	line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 2);
	hold off
	set(gca, 'XScale', 'log');
	xlim([data_MTF.fms(1) data_MTF.fms(end)])
	xticks([2 5 10 20 50 100 200 500])
	xlabel('Modulation Frequency')
	ylim(MTF_range)
	set(gca,'fontsize',14)
	grid on
	box on
	%ylabel('Avg. Rate (sp/s)')
	legend('Unmodulated', 'Location','southwest', 'fontsize', 16, 'EdgeColor','none')
	if ineuron == 1
		title('MTF: BE', 'FontSize', 18)
	end

	%% STRF

	% Plot
	% nexttile(t, 3)
	% imagesc(data_STRF.t, data_STRF.f./1000, data_STRF.H2ex_strf-data_STRF.H2in_strf, data_STRF.clims_strf);
	% set(gca,'Ydir','normal','XLim',data_STRF.tlims, 'YLim',[0 plot_range(2)]./1000)
	% colormap(redblue)
	% xlabel('Time (s)');
	% ylabel('Frequency (kHz)')
	% title('STRF')
	% set(gca,'fontsize',12)


	%% Synth Timbre
	num_ds = size(data_ST, 1);
	if num_ds == 4
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
	else
		data_colors = {'#82BB95', '#03882F', '#034E1C'};
	end

	% Sort
	spls = cellfun(@(p) p.spl, params_ST);
	[~, order] = sort(spls);

	%nexttile(t)
	h(5) = subplot(1, 6, 5);
	hold on
	label_ind = 1;
	for ind = 1:num_ds

		rate = data_ST{order(ind)}.rate;
		rate_std = data_ST{order(ind)}.rate_std;
		rlb = data_ST{order(ind)}.rlb;
		rub = data_ST{order(ind)}.rub;
		fpeaks = data_ST{order(ind)}.fpeaks;
		spl = params_ST{order(ind)}.spl+3;

		% Plot
		rates_sm = smooth_rates(rate, rlb, rub);
		errorbar(fpeaks, rate, rate_std/sqrt(30), 'linewidth', 2, 'color', data_colors{ind})
	end

	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
	label(label_ind) = {'Estimated CF'};
	yline(data_RM.spont, 'Color','k', 'LineWidth',2)
	label(label_ind+1) = {'Spont'};
	xlabel('Spectral Peak  Frequency (Hz)')
	%ylabel('Avg. rate (sp/s)')
	set(gca, 'Fontsize', 14, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	xlim(plot_range);
	grid on
	box on
	if ineuron == 1
		title('Synthetic Timbre', 'FontSize',18);
	end

	for ind = 1:num_ds
		annotation('textbox',[0.655 0.82-0.051*(ind-1) 0.0829 0.0869],...
			'Color',data_colors{ind},...
			'String',[num2str(params_ST{order(ind)}.spl) ' dB SPL'],...
			'FontSize',16,...
			'EdgeColor','none');
	end

	% Plot 1/2 height bandwidth
	% dsids = cellfun(@(p)p.dsid, params);
	% datasets = categorical(dsids);
	% nexttile(3)
	% bar(datasets, width)
	% title('Half-Height BW')
	% xlabel('Datasets')
	% ylabel('BW (Hz)')
	% grid on
	% box on
	% ylim([0 max(width+20)])
	% set(gca,'FontSize',5)

	%% Run SFIE Model

	% if ismac
	% 	addpath '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Modeling/Original SFIE'
	% 	addpath '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/UR_EAR_2022a/source'
	% else
	% 	addpath 'C:\Users\jfritzinger\Box\02 - Code\Modeling\Original SFIE'
	% end
	% 
	% % Set up stimulus
	% for ind = 1:num_ds
	% 	params_ST{ind}.Fs = 100000;
	% 	params_ST{ind}.physio = 1;
	% 	params_ST{ind}.mnrep = 20;
	% 	params_ST{ind}.dur = 0.3;
	% 	[params_ST{ind}] = generate_synthetictimbre(params_ST{ind});
	% 	params_ST{ind}.num_stim = size(params_ST{ind}.stim, 1);
	% end
	% 
	% % Run model
	% for ist = 1:num_ds
	% 	timerVal = tic;
	% 
	% 	% Model parameters
	% 	model_params.type = 'SFIE';
	% 	model_params.range = 2; % 1 = population model, 2 = single cell model
	% 	model_params.species = 1; % 1 = cat, 2 = human
	% 	model_params.BMF = 100;
	% 	model_params.CF_range = CF;
	% 	model_params.num_CFs = 1;
	% 	model_params.CFs = CF;
	% 	model_params.nAN_fibers_per_CF = 5;
	% 	model_params.cohc = 1; % (0-1 where 1 is normal)
	% 	model_params.cihc = 1; % (0-1 where 1 is normal)
	% 	model_params.nrep = 3; % how many times to run the AN model
	% 	model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
	% 	model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
	% 	model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
	% 	model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
	% 	model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
	% 
	% 	% Run model
	% 	AN_HSR = modelAN(params_ST{ist}, model_params); % HSR for IC input
	% 	SFIE{ist} = wrapperIC(AN_HSR.an_sout, params_ST{ist}, model_params); % SFIE output
	% 
	% 	elapsedTime = toc(timerVal)/60;
	% 	disp(['SFIE model took ' num2str(elapsedTime) ' minutes'])
	% end
	% 
	% filename = [session '_Model.mat'];
	% save(fullfile(path, filename), 'params_ST', 'SFIE', 'model_params', 'AN_HSR', '-v7.3')

	filename = [session '_Model.mat'];
	load(fullfile(base, fpath, filename))

	% Plot SFIE model
	if num_ds == 4
		model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};
	else
		model_colors = {'#E29A63', '#984102', '#4E2201'};
	end

	% Sort
	spls = cellfun(@(p) p.spl, params_ST);
	[~, order] = sort(spls);

	%nexttile(t)
	h(6) = subplot(1, 6, 6);
	hold on
	for ind = 1:num_ds
		[rate, ~] = plotST(params_ST{order(ind)}, SFIE{order(ind)}.average_ic_sout_BE, 0);
		plot(params_ST{order(ind)}.fpeaks,rate, ...
			'linewidth', 2, 'Color', model_colors{ind});
		R_int = corrcoef(data_ST{order(ind)}.rate,rate);
		R(ind) = R_int(1, 2).^2;
	end
	xlabel('Spectral Peak Frequency (Hz)')
	ylim(ylimits)
	set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	xlim(plot_range);
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
	set(gca,'FontSize',14)
	box on
	grid on
	%ylabel('Avg. Rate (sp/s)')
	if ineuron == 1
		title('SFIE Model Prediction', 'FontSize',18)
	end

	% R^2 annotation
	for ind = 1:num_ds
		if ineuron == 1
			annotation('textbox',[0.85 0.83-0.051*(ind-1) 0.0829 0.0917],...
				'Color',model_colors{ind},...
				'String',['R^2 = ' num2str(round(R(ind), 2))],...
				'FontSize',16,...
				'FitBoxToText','off',...
				'EdgeColor','none');
			annotation('textbox',[0.767 0.82-0.051*(ind-1) 0.0829 0.0869],...
				'Color',model_colors{ind},...
				'String',[num2str(params_ST{order(ind)}.spl) ' dB SPL'],...
				'FontSize',16,...
				'EdgeColor','none');
		else
			annotation('textbox',[0.85 0.23-0.051*(ind-1) 0.0829 0.0917],...
				'Color',model_colors{ind},...
				'String',['R^2 = ' num2str(round(R(ind), 2))],...
				'FontSize',16,...
				'FitBoxToText','off',...
				'EdgeColor','none');
			annotation('textbox',[0.767 0.22-0.051*(ind-1) 0.0829 0.0869],...
				'Color',model_colors{ind},...
				'String',[num2str(params_ST{order(ind)}.spl) ' dB SPL'],...
				'FontSize',16,...
				'EdgeColor','none');
		end
	end


	%% Rearrange

	% % Horizontal
	% set(h(1), 'position', [0.043,0.1312,0.1946,0.2604])
	% set(h(2), 'position', [0.043,0.3916,0.1946,0.2604])
	% set(h(3), 'position', [0.043,0.6520,0.1946,0.2604])
	% set(h(4), 'position', [0.282,0.1312,0.1946,0.7813])
	% set(h(5), 'position', [0.522,0.1312,0.1946,0.7813])
	% set(h(6), 'position', [0.762,0.1312,0.1946,0.7813])

	set(h(1), 'position', [0.043,0.1312,0.21,0.2604])
	set(h(2), 'position', [0.043,0.3916,0.21,0.2604])
	set(h(3), 'position', [0.043,0.6520,0.21,0.2604])
	set(h(4), 'position', [0.282,0.1312,0.21,0.7813])
	set(h(5), 'position', [0.522,0.1312,0.21,0.7813])
	set(h(6), 'position', [0.762,0.1312,0.21,0.7813])


	% RM
	annotation('textbox',[0.049 0.284 0.082 0.091],...
		'Color',[0.643137254901961 0.607843137254902 0.815686274509804],...
		'String','30 dB SPL',...
		'FontSize',16,...
		'FitBoxToText','off',...
		'EdgeColor','none');
	annotation('textbox',[0.0466 0.811 0.082 0.0917],...
		'Color',[0.125490196078431 0.0666666666666667 0.419607843137255],...
		'String','70 dB SPL',...
		'FontSize',16,...
		'FitBoxToText','off',...
		'EdgeColor','none');
	annotation('textbox',[0.047 0.550 0.082 0.091],...
		'Color',[0.368627450980392 0.313725490196078 0.662745098039216],...
		'String','50 dB SPL',...
		'FontSize',16,...
		'FitBoxToText','off',...
		'EdgeColor','none');
	if ineuron == 1
		annotation('textbox',[0.194 0.281 0.028 0.0917],...
			'Color',[0.301960784313725 0.301960784313725 0.301960784313725],...
			'String','CF',...
			'FontSize',16,...
			'FitBoxToText','off',...
			'EdgeColor','none');
	else
		annotation('textbox',[0.21 0.281 0.028 0.0917],...
			'Color',[0.301960784313725 0.301960784313725 0.301960784313725],...
			'String','CF',...
			'FontSize',16,...
			'FitBoxToText','off',...
			'EdgeColor','none');
	end

end

