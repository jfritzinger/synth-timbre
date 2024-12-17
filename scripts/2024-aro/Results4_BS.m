%% R024S470, 8/20/2021, TT2N1, Cell: Excellent, CF = 1170Hz, Depth 7.75
close all
clear

%% Phyiso session matched to model

rab_num = 24;
session = 'R024S470';
TT = 2;
N = 1;
CF = 1200;

%% Load in physio & analyze
% Identify current machine, which changes paths
[userid, base_dir, ~, report_path, data_path] = findPaths();
if ismac
	session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
else
	session_dir_name = base_dir{contains(base_dir, rab_num)};
end
session_dir = fullfile(session_dir_name, session);

% Load in physio session
[clusters, params, stims] = loadPhysiologySession(session_dir, session, userid);
cluster = clusters([clusters.tetrode] == TT); % Select tetrode
cluster = cluster([cluster.neuron] == N); % Select neuron

% Analysis
has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&...
	strcmp(p.SPEC_slide_type,'Spectral_Centroid')&&...
	p.fpeak_mid==1000,params);
spl = cellfun(@(p) p.spl,params(has_sc));
[~,sorted_spl] = sort(spl);

num_DSIDs = sum(has_sc);
data = cell(num_DSIDs, 1);
params = params(has_sc);
params = params(sorted_spl);
[params, fig, data] = plotPhysST(cluster, params, stims, CF, data, []);

%% Run models

% Set up stimulus
for ind = 1:4
	params{ind}.Fs = 100000;
	params{ind}.physio = 1;
	params{ind}.mnrep = 10;
	params{ind}.dur = 0.3;
	[params{ind}] = generate_synthetictimbre(params{ind});
	params{ind}.num_stim = size(params{ind}.stim, 1);
end

% Set up cell arrays for model results
num_MTF = 1;
num_stim = 4;
energy = cell(num_stim, num_MTF);
SFIE = cell(num_stim, num_MTF);
SFIE_eff = cell(num_stim, num_MTF);
lateral_model = cell(num_stim, num_MTF);
for iMTF = 1:num_MTF
	for ist = 1:num_stim

		%% Energy model
		timerVal = tic;

		CFs = CF;
		Fs = params{ist}.Fs;
		stimulus = [params{ist}.stim zeros(size(params{ist}.stim,1),0.1*Fs)];
		gamma_param.srate = Fs;
		tvals = (1:length(stimulus))/Fs;
		gamma_IF_reg = zeros(length(CFs),length(tvals));
		impaired = 0; % 0 = not impaired; 1 = 'impaired'

		pin_gamma = zeros(size(stimulus, 1), Fs*params{ist}.dur+0.1*Fs);
		for istim = 1:size(stimulus, 1)
			gamma_param.fc = CFs;
			pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
		end
		pin_gamma = pin_gamma(:,1:params{ist}.dur*Fs);
		energy{ist} = sqrt(mean(pin_gamma.^2,2));

		elapsedTime = toc(timerVal)/60;
		disp(['Energy model took ' num2str(elapsedTime) ' minutes'])

		%% Set up SFIE model and run
		timerVal = tic;

		% Model parameters
		model_params.type = 'SFIE';
		model_params.range = 2; % 1 = population model, 2 = single cell model
		model_params.species = 1; % 1 = cat, 2 = human
		model_params.BMF = 100;
		model_params.CF_range = CF;
		model_params.num_CFs = 1;
		model_params.CFs = CF;
		model_params.nAN_fibers_per_CF = 5;
		model_params.cohc = 1; % (0-1 where 1 is normal)
		model_params.cihc = 1; % (0-1 where 1 is normal)
		model_params.nrep = 3; % how many times to run the AN model
		model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
		model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
		model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
		model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
		model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

		% Run model
		AN_HSR = modelAN(params{ist}, model_params); % HSR for IC input
		SFIE{ist, iMTF} = wrapperIC(AN_HSR.an_sout, params{ist}, model_params); % SFIE output

		elapsedTime = toc(timerVal)/60;
		disp(['SFIE model took ' num2str(elapsedTime) ' minutes'])

		%% Set up SFIE efferent model and run
		timerVal = tic;
		clear model_params

		model_params.type = 'SFIE';
		model_params.range = 2; % 1 = population model, 2 = single cell model
		model_params.species = 1; % 1 = cat, 2 = human
		model_params.BMF = 100;
		model_params.CF_range = CF;
		model_params.num_CFs = 1;
		model_params.CFs = CF;
		model_params.nAN_fibers_per_CF = 10;
		model_params.cohc = 1; % (0-1 where 1 is normal)
		model_params.cihc = 1; % (0-1 where 1 is normal)
		model_params.nrep = 3; % how many times to run the AN model
		model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
		model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
		model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
		model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050

		nstim = size(params{ist}.stim,1);
		CFs = CF;
		nCFs = length(CFs);
		SFIE_eff{ist, iMTF}.average_ic_sout_BE = zeros(nstim, nCFs);
		AN_HSR.average_AN_sout = zeros(nstim, nCFs);
		AN_LSR.average_AN_sout = zeros(nstim, nCFs);
		stim_dur = size(params{ist}.stim, 2);
		ihc_eff = zeros(nstim, 1, stim_dur);
		hsr_eff = zeros(nstim, 1, stim_dur);
		lsr_eff = zeros(nstim, 1, stim_dur);
		ic_eff = zeros(nstim, 1, stim_dur);
		gain_eff = zeros(nstim, 1, stim_dur);
		for istim = 1:nstim
			stimulus = params{ist}.stim(istim,:);

			% Call model w/ efferent system enabled
			[ihc_eff(istim,:,:), hsr_eff(istim,:,:), lsr_eff(istim,:,:), ...
				ic_eff(istim,:,:), gain_eff(istim,:,:)] = sim_efferent_model( ...
				stimulus, CFs, moc_weight_wdr=2.0, moc_weight_ic=8.0, species=1);


			% Average and store IC rate
			AN_HSR.average_AN_sout(istim,:) = mean(hsr_eff(istim,:,:),3);
			AN_LSR.average_AN_sout(istim,:) = mean(lsr_eff(istim,:,:),3);
		end
		SFIE_eff{ist, iMTF} = wrapperIC(hsr_eff, params{ist}, model_params); % SFIE output

		elapsedTime = toc(timerVal)/60;
		disp(['SFIE with efferents model took ' num2str(elapsedTime) ' minutes'])

		%% Set up lateral inhibition model and run
		% timerVal = tic;
		% clear model_params
		%
		% % Lateral Model
		% model_params.type = 'Lateral Model';
		% model_params.lateral_CF = [CF/2 CF CF*2];
		% model_params.CFs = model_params.lateral_CF;
		% model_params.CF_range = model_params.CFs(2);
		% if iMTF == 1 % BE
		% 	CS_param = [0.35, 0.35, 0.003];
		% 	model_params.config_type = 'BS inhibited by off-CF BS';
		% else % BS
		% 	CS_param = [0.08, 0.08, 0.003];
		% 	model_params.config_type = 'BS inhibited by off-CF CN';
		% end
		%
		% % Model parameters
		% model_params.range = 2; % 1 = population model, 2 = single cell model
		% model_params.species = 1; % 1 = cat, 2 = human
		% model_params.num_CFs = 1;
		% model_params.nAN_fibers_per_CF = 10;
		% model_params.cohc = 1; % (0-1 where 1 is normal)
		% model_params.cihc = 1; % (0-1 where 1 is normal)
		% model_params.nrep = 1; % how many times to run the AN model
		% model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
		% model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
		% model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
		% model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
		% model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
		% model_params.BMF = 100;
		%
		% % Run Lateral Inhibition Model
		% AN = modelLateralAN(params{ist}, model_params);
		% lateral_model{ist, iMTF} = modelLateralSFIE(params{ist}, ...
		% 	model_params, AN.an_sout, AN.an_sout_lo, AN.an_sout_hi,'CS_params', CS_param);
		%
		% elapsedTime = toc(timerVal)/60;
		% disp(['Lateral inhibition model took ' num2str(elapsedTime) ' minutes'])
	end
end

%% Figure Properties

data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};

ylimits = [0 80];
figure('Position', [460,427,1300,420])
tiledlayout(1, 4, 'TileSpacing','compact')


% Plot data
clear leg
nexttile
hold on
for ind = 1:4
	errorbar(data{ind}.fpeaks,data{ind}.rate,data{ind}.rate_std/(sqrt(params{ind}.nrep)),...
		'LineWidth', 1, 'Color', data_colors{ind});
	leg(ind,:) = [num2str(params{ind}.spl+3) ' dB SPL'];
end
%plot(data{ind}.fpeaks{1}([1 end]),[1 1]*data{ind}.spont,'-','Color',[0.5 0.5 0.5], 'LineWidth', 2)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
legend(char(leg, 'CF'),...
	'location', 'northeast', 'FontSize',12)

plot_range = [data{ind}.fpeaks(1) data{ind}.fpeaks(end)];
xlabel('Spectral Peak Frequency (Hz)')
xlim(plot_range);
set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
ylabel('Avg. Rate (sp/s)')
ylim(ylimits)
set(gca,'FontSize',14)
box on
grid on
title('Single-Unit', 'FontSize',18);


% Plot Models
% Energy
clear leg
nexttile
hold on
for ind = 1:4
	[rate, ~] = plotSyntheticTimbre(params{ind}, energy{ind}, 0);
	plot(params{ind}.fpeaks,rate, ...
		'linewidth', 1.5, 'Color', model_colors{ind});
	R_int = corrcoef(data{ind}.rate,rate);
	R(ind) = R_int(1, 2).^2;
	leg(ind,:) = [num2str(params{ind}.spl+3) ' dB SPL, R^2 = ' num2str(round(R(ind), 3))];
end
xlabel('Spectral Peak Frequency (Hz)')
%ylim([0 0.003])
set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
xlim(plot_range);
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
legend(char(leg,'CF'),'location', 'northeast', 'FontSize',12)
set(gca,'FontSize',14)
box on
grid on
title('Energy', 'FontSize',18)
ylabel('Avg. Rate (sp/s)')

% SFIE
clear leg
nexttile
hold on
for ind = 1:4
	[rate, ~] = plotSyntheticTimbre(params{ind}, SFIE{ind, iMTF}.average_ic_sout_BS, 0);
	plot(params{ind}.fpeaks,rate, ...
		'linewidth', 1.5, 'Color', model_colors{ind});
	R_int = corrcoef(data{ind}.rate,rate);
	R(ind) = R_int(1, 2).^2;
	leg{ind,:} = [num2str(params{ind}.spl+3) ' dB SPL, R^2 = ' num2str(round(R(ind), 3))];
end
leg{ind+1, :} = 'CF';
xlabel('Spectral Peak Frequency (Hz)')
ylim(ylimits)
set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
xlim(plot_range);
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
legend(leg,'location', 'northeast', 'FontSize',12)

set(gca,'FontSize',14)
box on
grid on
title('SFIE', 'FontSize',18)
ylabel('Avg. Rate (sp/s)')

% SFIE with efferent
clear leg
nexttile
hold on
for ind = 1:4
	[rate, rate_std] = plotSyntheticTimbre(params{ind}, SFIE_eff{ind, iMTF}.average_ic_sout_BS, 0);
	plot(params{ind}.fpeaks,rate, ...
		'linewidth', 1.5, 'Color', model_colors{ind});
	R_int = corrcoef(data{ind}.rate,rate);
	R(ind) = R_int(1, 2).^2;
	leg(ind,:) = [num2str(params{ind}.spl+3) ' dB SPL, R^2 = ' num2str(round(R(ind), 3))];
end
xlabel('Spectral Peak Frequency (Hz)')
ylim(ylimits)
set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
xlim(plot_range);
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
legend(char(leg,'CF'),'location', 'northeast', 'FontSize',12)

set(gca,'FontSize',14)
box on
grid on
title('SFIE with efferent', 'FontSize',18)
ylabel('Avg. Rate (sp/s)')

%% Plot all models on top of each other

ylimits = [0 80];
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
figure;
set(gcf, 'Position', [460,427,1300,420])
set(gcf,'color','w');
tiledlayout(1, 4, 'TileSpacing','compact')

for ind = 1:4

	% Data
	nexttile
	hold on
	plot(data{ind}.fpeaks,data{ind}.rate,'LineWidth', 1.5);
	leg{1}  = 'Data';

	% Energy
	[e_rate, ~] = plotSyntheticTimbre(params{ind}, energy{ind}, 0);
	plot(params{ind}.fpeaks,e_rate*400, 'linewidth', 1.5);
	R_int = corrcoef(data{ind}.rate,e_rate);
	R(ind) = R_int(1, 2).^2;
	leg{2}  = ['Energy, R^2 = ' num2str(round(R(ind), 3))];

	% SFIE
	[SFIE_rate, ~] = plotSyntheticTimbre(params{ind}, SFIE{ind, iMTF}.average_ic_sout_BS, 0);
	plot(params{ind}.fpeaks,SFIE_rate, 'linewidth', 1.5);
	R_int = corrcoef(data{ind}.rate,SFIE_rate);
	R(ind) = R_int(1, 2).^2;
	leg{3}  = ['SFIE, R^2 = ' num2str(round(R(ind), 3))];

	% SFIE with efferent
	[SFIEe_rate, rate_std] = plotSyntheticTimbre(params{ind}, SFIE_eff{ind, iMTF}.average_ic_sout_BS, 0);
	plot(params{ind}.fpeaks,SFIEe_rate, 'linewidth', 1.5);
	R_int = corrcoef(data{ind}.rate,SFIEe_rate);
	R(ind) = R_int(1, 2).^2;
	leg{4} = ['Efferent, R^2 = ' num2str(round(R(ind), 3))];


	% Parameters
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
	legend(leg,'location', 'northeast', 'FontSize',12)
	plot_range = [data{ind}.fpeaks(1) data{ind}.fpeaks(end)];
	xlabel('Spectral Peak Frequency (Hz)')
	xlim(plot_range);
	set(gca, 'XTick', plot_range(1)+200:400:plot_range(2)-200);
	ylabel('Avg. Rate (sp/s)')
	ylim(ylimits)
	set(gca,'FontSize',14)
	box on
	grid on
	title([num2str(params{ind}.spl+3) ' dB SPL'])
end

%% Plot MTF
%
% has_mtf = cellfun(@(d)strcmp(d.param.type,'typMTFN'), data_full);
% data = data_full(has_mtf);
% ind = 1;
%
% figure;
% set(gcf, 'Position', [560,527,600,400])
% set(gcf,'color','w');
%
% % Calculate MTF
% dur = data{ind}.param.dur/1000; % stimulus duration in seconds.
% all_mod_depths = double([data{ind}.param.list.mdepth]).';
% all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
% [~,~,mdi] = unique(all_mod_de	pths);
% [fms,~,fmi] = unique(double([data{ind}.param.list.fm]).');
% if fms(1) == 0
%     fms(1) = 1.2;
% end
%
% num_mod_freqs = length(fms);
% num_depths = length(data{ind}.param.all_mdepths);
% map_size = [num_mod_freqs, num_depths];
%
% spike_rates = data{ind}.spike_info.num_spikes_delayed/...
%     (dur - data{ind}.spike_info.onsetWin/1000);
% [rate_MTF,rate_std,rlb,rub] = accumstats({fmi,mdi},spike_rates, map_size);
% rate_sm = zeros(size(rate_MTF));
% for j = 1:num_depths
%     rate_sm(:,j) = smooth_rates(rate_MTF(:,j),rlb(:,j),rub(:,j));
% end
%
% % Plot
% hold on
% line([1 fms(end)], [1 1]*rate_MTF(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
% errorbar(fms,rate_MTF,rate_std,'Color','k', 'Marker','.', 'MarkerSize',10, 'LineWidth', 2);
% plot(fms,rate_sm,'-b', 'LineWidth', 2)
%
% % Label the plots
% box on
% grid on
% xtick = [1 2 5 10 20 50 100 200 500];
% xlim(xtick([1 end]))
% xlabel('Modulation Freq (Hz)')
% ylabel('Spike Rate (sp/s)')
% set(gca,'XTick',xtick,'XScale', 'log')
% hold off
% set(gca,'fontsize',20)
% legend('No modulation', 'fontsize',19, 'location', 'northwest')
% title('Band-Suppressed', 'FontSize', 24);
%
%
% %% RM
%
% figure;
% set(gcf, 'Position', [560,527,800,400])
% set(gcf,'color','w');
% hold on
%
% % Plot RM
% plot(RM.freqs,RM.rates(:,5),'color', '#20116B','LineWidth',2) % 70 dB
% plot(RM.freqs,RM.rates(:,4),'color', '#5E50A9','LineWidth',2) % 50 dB 7B6DC1
% plot(RM.freqs,RM.rates(:,3),'color', '#A49BD0','LineWidth',2) % 30 dB
% plot(RM.freqs([1 end]),[1 1]*RM.spont,'-','LineWidth',2, 'Color',[0.7 0.7 0.7])
% xline(BS.CF, '--', 'Color', [0.4 0.4 0.4],'LineWidth',3);
% box on
% grid on
% hold off
% %ylim([0 RM.max_y+10])
% %set(gca,'XTick',[])
% %xticks(x_label)
% %yticks([0 100 200])
% %xticklabels([])
% set(gca,'fontsize',20)
% title('Response Map', 'FontSize', 24)
% set(gca,'XTick',[250 500 1000 2000 5000 10000])
% xlim([250 10000])
% legend('70 dB SPL', '50 dB SPL', '30 dB SPL', 'Spont. Rate', 'CF', 'fontsize',19)
% ylabel('Avg. Rate (sp/s)')
% xlabel('Tone Frequency (Hz)')
% set(gca, 'XScale', 'log');

