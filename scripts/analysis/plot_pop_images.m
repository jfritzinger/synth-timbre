%% plot_pop_images.m
%
% Script that plots all sythetic timbre responses as images aligned by CF.
% Neurons can be sorted by CF, MTF type, or peak/dip response. 
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

% Load in spreadsheet with peak information
spreadsheet_name = 'peak_picking.xlsx';
table = readtable(fullfile(datapath, spreadsheet_name));

% Find number of neurons in each category 
% Find sessions for target MTF type
%MTF_target = 'F';
%isMTF = strcmp(sessions.MTF, MTF_target);
%isMTF = contains(sessions.MTF, 'H');

% Find sessions for target synthetic timbre response
% bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
% bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
% bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
% bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
% bin200_MTF = bin200; % & isMTF;
%has_data = bin200_MTF(:,1) | bin200_MTF(:,2) | bin200_MTF(:,3) | bin200_MTF(:,4);
%indices = find(has_data);
%num_index = length(indices);

%% Plot imagesc of all BS responses sorted by CF

figure('position', [60,30,966,867])
backgroundcolor =  [0.8 0.8 0.8];
spl = [43, 63, 73, 83];
spls = {'43', '63', '73', '83'};


for ispl = 1:4

	% Find peaks and dips from table 
	isspl = table.SPL == spl(ispl);
	ispeak = strcmp(table.Type, 'Flat');
	is200 = table.F0 == 200;
	isbin = table.binmode == 2;
	isall = isspl &  ispeak & is200 & isbin;
	putatives = table.Putative(isall);
	num_index = size(putatives, 1);

	CFs = table.CF(isall);
	CF_names = cell(num_index, 1);
	array_z = NaN(num_index,10000);

	% %indices = find(bin200_MTF(:,ispl));
	% %num_index = length(indices);
	% array_z = NaN(num_index,10000);
	% %CFs = NaN(num_index, 1);
	% CFs = sessions.CF(indices);
	% CF_names = cell(num_index, 1);

	for isesh = 1:num_index
		%if bin200_MTF(indices(isesh), ispl)==1

			% Load in session
			%putative = sessions.Putative_Units{indices(isesh)};
			%CF = sessions.CF(indices(isesh));
			putative = putatives{isesh};
			CF = CFs(isesh);
			load(fullfile(datapath, 'neural_data', [putative '.mat']))
			params_ST = data(5+ispl, 2);

			% peak_ind = find(strcmp(putative, table.Putative) & spl(ispl)==table.SPL);
			% if ~isempty(peak_ind)
			% 	type = table.Type{peak_ind};
			% 	if strcmp(type, 'Peak')||strcmp(type, 'Dip')
				% 	continue
			% 	end
			% end

			%CFs(isesh) = CF;
			CF_names{isesh} = [num2str(round(CFs(isesh))) ' Hz'];

			% Analysis
			data_ST = analyzeST(params_ST, CF);
			data_ST = data_ST{1};
			params_RM = data{2,2};
			data_RM = analyzeRM(params_RM);			
			spont = data_RM.spont;

			% General analysis
			rate = data_ST.rate;
			rate = rate - spont;
			fpeaks = data_ST.fpeaks;
			fpeaks_re_CF = log2(fpeaks/CF);
			%fpeaks_re_CF = fpeaks/CF;

			% Align by CF (approximately)
			f = linspace(-3, 3, 10000);
			[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
			[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
			f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

			% Interpolate & get z-score
			r_interp = interp1(fpeaks_re_CF, rate,f_interp, 'spline');
			z_rate = zscore(r_interp);
			%z_rate = r_interp;
			array_z(isesh, f_ind(1):f_ind(2)-1) = z_rate;
		%else
		%	CF_names{isesh} = [num2str(round(CFs(isesh))) ' Hz'];
		%end
	end

	% Order by CF
	[~, max_ind] = sort(CFs);
	CF_order = array_z(max_ind,:);
	CFs_ordered = CF_names(max_ind);

	%nexttile
	s(ispl) = subplot(1, 4, ispl);
	h = imagesc(f, 1:size(CF_order, 1), CF_order);
	xline(0, 'k')
	%colormap(redblue)
	set(h, 'AlphaData', ~isnan(CF_order))
	set(gca,'color',backgroundcolor);
	if ispl == 1
		ylabel('Neuron Number', 'Color','w')
	end
	if ispl == 2
		xlabel('Spectral Peak Freq w.r.t. CF (Oct.)')
	end
	yticklabels([])
	xlim([-2 2])
	clim([-2 2])
	xticklabels({'-2', '-1', 'CF', '1', '2'})
	set(gca, 'fontsize', 14)
	title([spls{ispl} ' dB SPL'], 'fontsize', 18)
end

for ind = 1:size(CF_order, 1)
	CF_label = CFs_ordered{ind};
	interval = 0.91 / size(CF_order, 1);
	start = 0.01;
	annotation('textbox',[start 0.935-interval*ind 0.07 0.049], ...
		'String',CF_label, 'EdgeColor',...
		'none', 'FontSize',8);
end

set(s(1), 'position', [0.06,0.0552, 0.205,0.9071])
set(s(2), 'position', [0.29,0.0552,0.205,0.9071])
set(s(3), 'position', [0.52,0.0552,0.205,0.9071])
set(s(4), 'position', [0.75,0.0552,0.245,0.9071])
colorbar
