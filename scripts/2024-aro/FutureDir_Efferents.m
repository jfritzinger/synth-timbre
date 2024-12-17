%% FutureDir_Efferents
% J. Fritzinger, 1/16/24
clear

%% Physio Example

figure('Position',[343,611,850,259])
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact')

nexttile
hold on

% R027S058, TT3N2
% session = 'R027S058_TT3_N2';
% CF = 2297;
% ds = 9;

% R024S473, 
% session = 'R024S473_TT2_N1';
% CF = 1150;
% ds = 11;

% R025S688
session = 'R025S688_TT4_N1';
CF = 2254;
ds = 10;

% Load in examples
base = getPaths();
fpath = 'data/aro-2024';
filename = fullfile(base, fpath, [session '.mat']);
load(filename, 'params', 'data', 'cluster', 'stims')

% Synth Timbre
data_ST = data{ds}; %
param = params{ds};

% Analysis
num_bins = 3;
start_bin = 0;
end_bin = 200;
win_start = linspace(start_bin, end_bin, num_bins);
win_end = win_start + 100;
color_gradient = {'#648FFF', '#FFB000', '#DC267F'};
num_win = length(win_start);
stimx = stims(ds);
this_ds = stims(ds).dsid == ds;
dur = param.dur/1000; % stimulus duration in seconds.

% Check if this session has an spl shift
spl_shift = checkIfSPLShift(param.animal, param.session);
spl = param.spl + spl_shift;

[fpeaks,~,fi] = unique([param.list.fpeak].');
num_fpeaks = length(fpeaks);
stim_time_bins = [stimx.times;[-1 2]*stimx.times(end-1:end)];
[~,abs_stim_num] = histc(cluster.t_spike,stim_time_bins);
valid = abs_stim_num ~= 0;
t_spike_rel = zeros(size(cluster.t_spike));
t_spike_rel(valid) = cluster.t_spike(valid) - stimx.times(abs_stim_num(valid));
rel_id = abs_stim_num;
rep = cell2mat({param.list.rep})';

% Sort
width = zeros(num_win, 1);
max_rate = zeros(num_win, 1);
for iwin = 1:num_win
	avg_rate_win = zeros(num_fpeaks, 1);
	std_rate_win = zeros(num_fpeaks, 1);
	for j2 = 1:num_fpeaks
		k = find(fi == j2);
		x = t_spike_rel(ismember(rel_id,k));
		y = rep(rel_id(ismember(rel_id,k)));
		x = x/1000; % ms
		win_y = y(x>win_start(iwin) & x<win_end(iwin));
		rate_win = histcounts(win_y,param.nrep);
		avg_rate_win(j2) = mean(rate_win)/0.05;
		std_rate_win(j2) = std(rate_win/0.05);
		lb(j2) = min(rate_win)/0.05;
		ub(j2) = max(rate_win)/0.05;

	end
	max_rate(iwin) = round(max(avg_rate_win, [],'all'))+5;

	% Smooth
	rates_sm = smooth_rates(avg_rate_win, lb, ub);
	max_rate(iwin) = round(max(rates_sm, [],'all'))+3;

	% Calculate half height width
	width(iwin) = findHalfHeightWidth(fpeaks, rates_sm, 0);

	% Plot
	rate_win_all(iwin, :) = avg_rate_win;
	plot(fpeaks/1000,avg_rate_win, 'LineWidth',2, 'Color',color_gradient{iwin});
end

xline(CF/1000, '--', 'Color', [0.5 0.5 0.5], 'LineWidth',2);
for ii = 1:num_win
	leg(ii,:) = {sprintf('%d-%dms',win_start(ii), win_end(ii))};
end

xlabel('Frequency (kHz)')
ylabel('Spike rate (sp/s)')
set(gca,'FontSize',16)
%legend(leg, 'Location', 'northwest')
grid on
title('BS Neuron, Time Lapse', 'fontsize', 18);
box on
plot_range = [param.fpeaks(1) param.fpeaks(end)];
xlim(plot_range/1000)

window_names = {'0-100 ms', '100-200 ms', '200-300 ms'};
for ind = 1:3 
	annotation('textbox',[0.09 0.78-0.07*(ind-1) 0.15 0.1],...
		'Color',color_gradient{ind},...
		'String',window_names{ind},...
		'FontSize',16,...
		'EdgeColor','none');
end

%% Model Example
% 
% %compile_mex
% 
% % Stimuli
% param.Fs = 100000;
% param.physio = 1;
% param.mnrep = 30;
% param.dur = 0.3;
% %param.spl = 70;
% [param] = generate_synthetictimbre(param);
% param.num_stim = size(param.stim, 1);
% 
% % Model parameters
% model_params.type = 'SFIE';
% model_params.range = 2; % 1 = population model, 2 = single cell model
% model_params.species = 1; % 1 = cat, 2 = human
% model_params.BMF = 100;
% CF_lo = CF * 2^(-3/4);
% CF_hi = CF * 2^(3/4);
% model_params.CF = CF;
% model_params.CF_range = CF;
% %model_params.CF_range = logspace(log10(CF_lo), log10(CF_hi), 21);
% model_params.num_CFs = 1;
% model_params.CFs = CF;
% odel_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
% model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
% 
% %nCF = 21;
% nCF = 1;
% nstim = size(param.stim,1);
% stim_dur = size(param.stim, 2);
% hsr_eff = zeros(nstim, nCF, stim_dur);
% ic_eff = zeros(nstim, nCF, stim_dur);
% gain_eff = zeros(nstim, nCF, stim_dur);
% %CFs = logspace(log10(CF_lo), log10(CF_hi), 21);
% CFs = CF;
% for istim = 1:nstim
% 	stimulus = param.stim(istim,:);
% 
% 	% Call model w/ efferent system enabled
% 	[~, hsr_eff(istim,:,:), ~, ic_eff(istim,:,:), ...
% 		gain_eff(istim,:,:)] = sim_efferent_model( ...
% 		stimulus, CFs, species=1, moc_weight_wdr=2.0,...
% 		moc_weight_ic=8.0);
% 
% end
% %SFIE_eff = wrapperIC(hsr_eff(:,11,:), param, model_params); % SFIE output
% SFIE_eff = wrapperIC(hsr_eff, param, model_params); % SFIE output
% 
% 
% save(fullfile(modelpath, 'FutureDir_Efferents.mat'), 'SFIE_eff', 'param')

%% Load model results 

load(fullfile(base, fpath, 'FutureDir_Efferents.mat'), 'SFIE_eff', 'param')

%% Plot
ic_BS = squeeze(SFIE_eff.ic_BS);
ic_BE = squeeze(SFIE_eff.ic_BE);

[fpeaks,~,fpeaksi] = unique([param.mlist.fpeak].');
num_fpeaks = length(fpeaks);

nexttile
hold on
for ii = 1:3 
	window = (ii-1)*param.Fs/10+1:ii*param.Fs/10;
	rates = ic_BS(:,window);
	avg_rates = mean(rates, 2);

	rate_size = [num_fpeaks,1];
	[rate,rate_std] = accumstats({fpeaksi},avg_rates, rate_size);
	plot(fpeaks/1000,rate,'color', color_gradient{ii}, 'LineWidth',2);

	R_int = corrcoef(rate_win_all(ii, :),rate);
	R(ii) = R_int(1, 2).^2;

end
xline(CF/1000, '--', 'color', [0.5 0.5 0.5], 'LineWidth',2)
xlabel('Frequency (kHz)')
%ylabel('Spike rate (sp/s)')
set(gca,'FontSize',16)
%legend(leg, 'Location', 'northwest')
grid on
title('BS Efferent Model', 'fontsize', 18);
box on
xlim(plot_range/1000)

% Annotation
window_names = {'0-100 ms', '100-200 ms', '200-300 ms'};
for ind = 1:3 
	annotation('textbox',[0.56 0.78-0.07*(ind-1) 0.3 0.1],...
		'Color',color_gradient{ind},...
		'String',[window_names{ind} ', R^2 = ' num2str(round(R(ind), 2))],...
		'FontSize',16,...
		'EdgeColor','none');
end

% nexttile 
% gain = squeeze(gain_eff);
% imagesc(gain(1:40, :))
% colorbar
% title('Gain', 'fontsize', 18)
% xlabel('Time')
% set(gca,'FontSize',16)
% ylabel('Spectral Peak #')

