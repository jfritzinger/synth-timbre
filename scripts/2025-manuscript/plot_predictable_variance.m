%% plot_predictable_variance.m
%
% Script that plots predictable variance calculation for all binaural
% responses to synthetic timbre at the four different levels. 
%
% Author: J. Fritzinger
% Created: -------; Last revision: 2024-10-22
%
% -------------------------------------------------------------------------
clear

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% % Find number of neurons in each category 
% % Find sessions for target MTF type
% MTF_target = 'BS';
% isMTF = strcmp(sessions.MTF, MTF_target);
% 
% % Find sessions for target synthetic timbre response
% bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
% bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
% bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
% bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
% bin200_MTF = bin200 & isMTF;

%% Only get sessions with synthetic timbre 

synth(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
synth(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
synth(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
synth(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);

synth(:,5) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_100);
synth(:,6) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_100);
synth(:,8) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_100);

synth(:,9) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
synth(:,10) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
synth(:,11) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
synth(:,12) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);

any_synth = any(synth(:,1:12), 2);
table = sessions(any_synth, :);

%% Calculate predictable variance for all neurons 

num_index = size(table,1);
V_p = NaN(num_index, 4);
for isesh = 1:num_index

	% Load in session
	putative = table.Putative_Units{isesh};
	CF = table.CF(isesh);
	load(fullfile(datapath, [putative '.mat']))

	for ispl = 1:4
		% Analysis
		params_ST = data(5+ispl, 2);
		if ~isempty(params_ST{1})
			data_ST = analyzeST(params_ST);
			data_ST = data_ST{1};
			V_p(isesh, ispl) = data_ST.V_p;
		end
	end
end

%% Plot
figure('position', [519,533,1150,305])
tiledlayout(1, 4, 'TileSpacing','tight', 'Padding','compact')
spls = [43, 63, 73, 83];
for ispl = 1:4
	
	nexttile(ispl)
	histogram(V_p(:,ispl), 20)
	%scatter(CFs/1000, V_p(:,ispl), 'filled')
	number = V_p(:,ispl);
	number(isnan(number)) = [];
	title([num2str(spls(ispl)) ' dB SPL, n=' num2str(length(number))])
	if ispl == 1
		ylabel('# neurons')
	end
	set(gca, 'fontsize', 17)
	grid on
	xlabel('Predictable Variance')
	xlim([0 1])
	ylim([0 100])
end

%% Export 

exportgraphics(gcf, fullfile(savepath, 'manuscript', 'predictable-variance.png'), 'Resolution', 600)