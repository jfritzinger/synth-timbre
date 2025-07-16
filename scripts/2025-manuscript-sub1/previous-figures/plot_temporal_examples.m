%% plot_temporal_examples

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


%% Load in example

figure('Position',[67,297,8.5*ppi,6*ppi])
tiledlayout(4, 5, "TileSpacing","compact", 'TileIndexing','columnmajor')
legsize = 9;
fontsize = 11;
labelsize = 18;

%%
ind = [1, 6, 11, 16, 2, 7, 12, 17, 3, 8, 13, 18, 4, 9, 14, 19, 5, 10, 15, 20];
for ii = 1:20

	switch ii
		case 1 
			putative = 'R25_TT1_P8_N02'; % blurred at CF
		case 2
			putative = 'R25_TT1_P8_N01'; % blurred at CF
		case 3 
			putative = 'R25_TT4_P7_N01'; % blurred at CF
		case 4
			putative = 'R29_TT4_P2_N03'; % blurred at CF
		case 5
			putative = 'R29_TT4_P2_N09'; % 2 peaks 	
		case 6
			putative = 'R29_TT1_P3_N02'; % 2 peaks 
		case 7
			putative = 'R29_TT4_P2_N10'; % 2 peaks 
		case 8
			%putative = 'R25_TT2_P9_N07'; % 2 peaks 
			%putative = 'R29_TT1_P3_N01';
			%putative = 'R29_TT1_P2_N03';
			putative = 'R29_TT1_P3_N05';
		case 9
			putative = 'R29_TT1_P2_N03'; % 2 peaks, odd
		case 10 
			putative = 'R25_TT1_P8_N15'; %'R25_TT1_P8_N03'; % 2 peaks, odd 
		case 11
			putative = 'R25_TT1_P8_N04'; %'R29_TT4_P2_N16'; % 2 peaks, odd 
		case 12 
			putative = 'R25_TT2_P9_N02'; % % 2 peaks, odd
		case 13
			putative = 'R29_TT3_P2_N05'; % multiple
		case 14
			putative = 'R29_TT4_P3_N04'; % multiple
		case 15 
			putative = 'R29_TT3_P2_N04'; % multiple
		case 16
			putative = 'R25_TT3_P8_N06'; % multiple
		case 17
			putative = 'R29_TT4_P5_N02'; %'R24_TT1_P12_N01';
		case 18
			putative = 'R29_TT3_P5_N03'; % rate only
		case 19
			putative = 'R27_TT4_P8_N10';
		case 20
			putative = 'R27_TT2_P8_N05';
	end
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);
	MTF_shape = sessions.MTF{index};

	%% Analysis

	% Synthetic timbre analysis
	params = data(7, 2);
	params = params(~cellfun(@isempty, params));
	data_ST  = analyzeST(params, CF);
	data_ST = data_ST{1};
	param = params{1};
	temporal = analyzeST_Temporal(param, data_ST);

	% Plot as heatmap
	num_fpeaks = length(data_ST.fpeaks);
	max_rate = max(temporal.p_hist, [], 'all');
	freq_values = round(param.fpeaks);
	p_hist = temporal.p_hist;
	edges = temporal.t_hist;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	t = linspace(0, 5, size(p_hist,2));

	grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
	grayMap = flipud(grayMap);

	h(ind(ii)) = subplot(4, 5, ind(ii));
	hh = pcolor(t, data_ST.fpeaks./1000, p_hist);
	set(hh, 'EdgeColor', 'none');
	hold on
	yline(CF/1000, 'r', 'LineWidth',2)
	%colorbar;
	%axis square;
	colormap(grayMap);
	xlim([0 5])
	max_rate = max(p_hist, [], "all");
	clim([0 max_rate-max_rate*0.3])
	if ismember(ii, 3)
		ylabel('                                Spectral Peak Freq. (Hz)')
	end
	xticks(0:5)
	if ismember(ii, [4, 8, 12, 16, 20])
		xlabel('Period (ms)')
	else
		xticklabels([])
	end
	box off

	set(gca, 'fontsize', fontsize)
	if ii == 17
		hleg = legend('', 'CF', 'Location','northwest', 'fontsize', ...
			legsize, 'box', 'off');
		hleg.ItemTokenSize = [8, 8];
	end

	% if ii == 1
	% 	title('Blurred near CF')
	% elseif ii == 5
	% 	title('Two+ Peaks')
	% elseif ii == 9
	% 	title('Not uniform across freq')
	% elseif ii == 13
	% 	title('3-4 Peaks')
	% elseif ii == 17
	% 	title('No locking')
	% end

	msg = sprintf('%d sp/s', max_rate);
	text(0.58, 1, msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize, 'Color','r')


end

%% Set positions 

height = 0.19;
width = 0.14;
bottom = fliplr(linspace(0.1, 0.76, 4));
left = linspace(0.07, 0.85, 5);

left = repmat(left, 1, 4);
bottom = reshape(repmat(bottom, 5, 1), 1, []);

for ii = 1:20
	set(h(ii), 'Position', [left(ii) bottom(ii) width height])
end

% Annotations
left = linspace(0.02, 0.81, 5);
annotation('textbox',[left(1) 0.965 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) 0.965 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3) 0.965 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(4) 0.965 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(5) 0.965 0.0826 0.0385],'String',{'E'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');



