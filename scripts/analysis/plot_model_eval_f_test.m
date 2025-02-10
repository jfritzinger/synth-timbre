%% plot_model_eval_f_test

%% save_r2_table.m
%
% Script to create an excel spreadsheet with all R2 model results for each
% neuron at each level. Only for binaural stimuli. Used for
% 'plot_model_evaluations.m'. 
%
%
% Author: J. Fritzinger
% Created: 2022-09-13; Last revision: 2024-09-26 
%
% -------------------------------------------------------------------------
clear 

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%%  Synthetic Timbre Data Table 
 
% Initialize spreadsheet columns
modelpath = '/Volumes/Synth-Timbre/data/manuscript';

% Find sessions for target synthetic timbre response
bin200 = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
isMTF = strcmp(sessions.MTF, 'BE')|strcmp(sessions.MTF, 'BS');
bin200_MTF = bin200 & isMTF;

% Add R^2 values to the spreadsheet 
spls = [43, 63, 73, 83];
has_data = bin200_MTF;
indices = find(has_data);
num_index = length(indices);

p_s_e = NaN(num_index, 1);
p_l_e = NaN(num_index, 1);
p_l_s = NaN(num_index, 1);
p_s_pop = NaN(num_index, 1);
CF_all = NaN(num_index, 1);
ispl = 2;
MTFs = sessions.MTF(indices);
MTF_ind(1, :) = cellfun(@(m) strcmp(m, 'BE'), MTFs);
MTF_ind(2, :) = cellfun(@(m) strcmp(m, 'BS'), MTFs);
for isesh = 1:num_index

	% Load in data
	putative = sessions.Putative_Units{indices(isesh)};
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'SFIE_pop_model', [putative '_SFIE_pop.mat']), 'SFIE_pop')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Analysis
	sfie_temp = SFIE{ispl}.rate;
	energy_temp =  energy{ispl}.rate;
	lat_inh_temp =  lat_inh{ispl}.rate;
	sfie_pop_temp = SFIE_pop{ispl}.rate;
	CF_all(isesh) =  sessions.CF(indices(isesh));
	params_ST = data(7, 2);
	data_ST = analyzeST(params_ST, CF_all(isesh));
	data_ST = data_ST{1};

	% Normalize to 1 (will not use this when models are fit to rate correctly)
	sfie_temp = sfie_temp ./ max(sfie_temp);
	energy_temp = energy_temp ./ max(energy_temp);
	lat_inh_temp = lat_inh_temp ./ max(lat_inh_temp);
	sfie_pop_temp = sfie_pop_temp ./ max(sfie_pop_temp);
	rate = data_ST.rate ./ max(data_ST.rate);

	% figure;
	% hold on
	% plot(sfie_temp)
	% plot(energy_temp)
	% plot(lat_inh_temp)
	% plot(sfie_pop_temp)
	% plot(rate)

	% F-test SFIE/energy
	% [~,p,~,~] = vartest2(energy_temp,sfie_temp);
	p = model_ttest(rate, energy_temp, sfie_temp);
	p_s_e(isesh) = log(p);

	% F-test LatInh/energy
	%[~,p,~,~] = vartest2(energy_temp,lat_inh_temp);
	p = model_ttest(rate, energy_temp, lat_inh_temp);
	p_l_e(isesh) = log(p);
	
	% F-test LatInh/SFIE
	%[~,p,~,~] = vartest2(lat_inh_temp,sfie_temp);
	p = model_ttest(rate, lat_inh_temp, sfie_temp);
	p_l_s(isesh) = log(p);

	% F-test SFIE/population SFIE 
	%[~,p,~,~] = vartest2(sfie_pop_temp,sfie_temp);
	p = model_ttest(rate, sfie_pop_temp', sfie_temp);
	p_s_pop(isesh) = log(p);

	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end


%% Plot results 

figure('Position',[1796,680,1292,321])
tiledlayout(1, 4)

name = 'SFIE / Energy';
plotFtest(name, CF_all,MTF_ind, p_s_e)


name = 'Lat Inh / Energy';
plotFtest(name, CF_all, MTF_ind, p_l_e)

name = 'Lat Inh / SFIE';
plotFtest(name, CF_all, MTF_ind, p_l_s)

name = 'SFIE / SFIE population';
plotFtest(name, CF_all, MTF_ind, p_s_pop)


%% FUNCTIONS --------------------------------------------------------------

function plotFtest(name, CFs, iMTFs, vals)

nexttile
hold on
scatter(CFs(iMTFs(2,:))./1000, vals(iMTFs(2,:)), 'filled', ...
	"MarkerEdgeColor",'k', "MarkerFaceAlpha",0.5)
scatter(CFs(iMTFs(1,:))./1000, vals(iMTFs(1,:)), 'filled', ...
	"MarkerEdgeColor",'k', "MarkerFaceAlpha",0.5)
yline(log(0.05), '--')
set(gca, 'XScale', 'log')
ylabel('Log p-value from t-test of MSE')
xlabel('CF (kHz)')
set(gca, 'fontsize', 16)
title(name)
xticks([0.1 0.2 0.5 1 2 5 10])
grid on

% annotation
sig = sum(vals<log(0.05));
total = length(vals);
msg = sprintf('%d/%d units significant', sig, total);

text(0.05, 0.1, msg, 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)
end

% -------------------------------------------------------------------------

