%% Results 2: BE Population Analysis
% J. Fritzinger, updated 1/9/23
clear

%% Load in data
% 
% % BE
% %datafile = '2024-01-10_093623_Data_BE.mat';
% 
% % BS
% datafile = '2024-01-10_111039_Data_BS';
% 
% % Reads in spreadsheet & load data
% if ismac
% 	%filename = '/Volumes/CarneyLab/Rabbit_data/session_table.xlsx';
% 	filename = '/Volumes/Rabbit_data/session_table.xlsx';
% 	load(fullfile('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - Timbre/Data/', datafile))
% 
% else
% 	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
% 	load(fullfile('C:\Users\jfritzinger\Box\02 - Code\Physiology Data\', datafile))
% 
% end
% sessions = readtable(filename, 'PreserveVariableNames',true);
% 
% % Get sessions for BE
% population = population(~cellfun(@isempty, population));
% num_sessions = length(population);

%% Sessions of interest 


spreadsheet_name = 'SynthSessions.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

base = getPaths();
fpath = 'data/aro-2024';


% isBS = strcmp(sessions.MTF, 'BS');
% num_sessions = sum(isBS);
% BS_ind = find(isBS);

%% Find number of neurons in each category 

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BS');
num_sessions = sum(isBE);
BE_ind = find(isBE);

bin200_43 = 0;
bin200_63 = 0;
bin200_73 = 0;
bin200_83 = 0;
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(base, fpath, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest1 = find(bin & F0_200 & level_43);
	if any(interest1)
		bin200_43 = bin200_43+1;
	end

	interest2 = find(bin & F0_200 & level_63);
	if any(interest2)
		bin200_63 = bin200_63+1;
	end

	interest3 = find(bin & F0_200 & level_73);
	if any(interest3)
		bin200_73 = bin200_73+1;
	end

	interest4 = find(bin & F0_200 & level_83);
	if any(interest4)
		bin200_83 = bin200_83+1;
	end
end

%% Plot each dataset 
% 
% for isesh = 1:num_sessions
% in
% 	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
% 		= finddata(population{isesh,1}(:,1));
% 
% 	% Find dataset of interest
% 	interest43 = find(bin & F0_100 & level_43);
% 	interest63 = find(bin & F0_100 & level_63);
% 	interest73 = find(bin & F0_100 & level_73);
% 	interest83 = find(bin & F0_100 & level_83);
% 	interest = [interest43', interest63', interest73', interest83'];
% 	if any(interest)
% 		figure('Position',[292,274,1264,420])
% 		tiledlayout(1, 3)
% 
% 		x_label = [1000 2000 3000 4000 6000 8000];
% 
% 		% Plot RM
% 		data_RM = population{isesh,1}{4, 2};
% 		if ~isempty(data_RM)
% 			CF = data_RM.CF;
% 			nexttile(1)
% 			hold on
% 			spont_color = [0.4 0.4 0.4];
% 			CF_color = [0.3 0.3 0.3];
% 			plot(data_RM.freqs,data_RM.rates(:,5),'color', '#20116B','LineWidth',2) % 73 dB
% 			plot(data_RM.freqs,data_RM.rates(:,4),'color', '#5E50A9','LineWidth',2) % 53 dB
% 			plot(data_RM.freqs,data_RM.rates(:,3),'color', '#A49BD0','LineWidth',2) % 33 dB
% 			plot(data_RM.freqs([1 end]),[1 1]*data_RM.spont,'-','LineWidth',2, 'Color',spont_color)
% 			xline(CF, '--', 'Color', CF_color,'LineWidth',2);
% 			box on
% 			grid on
% 			hold off
% 			ylim([0 max(data_RM.rates, [], 'all')+10])
% 			set(gca,'XTick',[])
% 			xlim([250 14000])
% 			xticks(x_label)
% 			set(gca, 'XScale', 'log');
% 			set(gcf, 'color', 'w')
% 			set(gca,'fontsize',12)
% 			ylabel('Avg. Rate (sp/s)')
% 			legend('73dB SPL', '53dB SPL', '33dB SPL', 'Spont. Rate', 'fontsize',12, 'location', 'northwest')
% 			title('Response Map')
% 			xlabel('Frequency')
% 		end
% 
% 		% Plot MTF
% 		params_MTF = population{isesh,1}{5, 1};
% 		data_MTF = population{isesh,1}{5, 2};
% 		nexttile(2)
% 		hold on
% 		line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',2);
% 		errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',2, 'Color','k');
% 		line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 2);
% 		hold off
% 		set(gca, 'XScale', 'log');
% 		xlim([data_MTF.fms(1) data_MTF.fms(end)])
% 		xticks([1 2 5 10 20 50 100 200 500])
% 		xlabel('Modulation Frequency')
% 		%ylim(MTF_range)
% 		set(gca,'fontsize',12)
% 		grid on
% 		box on
% 		title(['MTF: ' name], 'FontSize', 12)
% 		ylabel('Avg. Rate (sp/s)')
% 		legend('Unmodulated', 'Location','southwest')
% 
% 
% 		num_interest = length(interest);
% 		for ind = 1:num_interest
% 			ids = interest(ind);
% 			data = population{isesh,1}{ids, 2};
% 			param = population{isesh,1}{ids, 1};
% 
% 			% Get CF rate
% 			[~, CF_ind] = min(abs(CF-data.fpeaks));
% 			CF_rate = data.rates_sm(CF_ind);
% 
% 			% Get spl 
% 			spl_shift = checkIfSPLShift(param.animal, param.session);
% 			spl = param.spl + spl_shift;
% 			leg{ind} = [num2str(spl) ' dB SPL'];
% 
% 			% Plot BE Synth Timbre
% 			nexttile(3)
% 			hold on
% 			plot(data.fpeaks,data.rates_sm, 'linewidth', 1.5);			
% 
% 		end
% 		leg{ind+1} = 'CF';
% 		legend(leg, 'Location','best')
% 		xline(CF, '--', 'Color',CF_color, 'linewidth', 1.5)
% 		xlabel('Frequency (Hz)')
% 		ylabel('Avg. Rate (sp/s)')
% 		box on
% 		title('Synth Timbre: F0=200, Bin')
% 		grid on
% 		%set(gca, 'fontsize', 14)
% 		clear leg
% 	end
% end

%% Seperate data into groups

figure('Position',[242,552,1454,366])
h1 = subplot(1, 4, 1:3);
hold on
count = 0;
ind = 1;
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	CF = sessions.CF(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = find(bin & F0_200 & level_83);
	if any(interest)
		data = data{interest};

		% Get CF rate
		[~, CF_ind] = min(abs(CF-data.fpeaks));
		CF_rate = data.rate(CF_ind);
		CF_sm = data.rates_sm(CF_ind);

		% Calculate HHBW
		width(ind) = data.width;
		ind = ind + 1;

		% Plot
		% nexttile([1 3])
		% hold on
		% plot(data.fpeaks,data.rate, 'linewidth', 1.5);
		% scatter(CF, CF_rate, 'filled', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')

		
		plot(h1, data.fpeaks,data.rates_sm, 'linewidth', 1.5);
		scatter(h1, CF, CF_sm, 'filled', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')

		count = count + 1;
	end
end

% Plotting Params
xlabel(h1, 'Frequency (Hz)')
ylabel(h1, 'Avg. Rate (sp/s)')
box on
title(h1, 'Synth Timbre Smoothed: BE, 83 dB SPL, F0=200, Bin')
set(h1, 'fontsize', 14)

% nexttile(2)
% xlabel('Frequency (Hz)')
% ylabel('Avg. Rate (sp/s)')
% box on
% title('Synth Timbre: BE, 73 dB SPL, F0=200, Bin')
% set(gca, 'fontsize', 14)

% Width figure
h2 = subplot(1, 4, 4);
histogram(h2, width, 10)
xlabel(h2, 'Half-Height BW (Hz)')
ylabel(h2, '# Neurons')
%title(h2, 'Synth Timbre: BE, 73 dB SPL, F0=200, Bin')
set(h2, 'fontsize', 14)

%% Overlayed 

figure('Position',[242,552,1246,366])
tiledlayout(1, 3, 'TileSpacing','compact')
for ispl = 1:3

	isBE = strcmp(sessions.MTF, 'BS');
	if ispl == 1
		isSPL = sessions.("43dB")==1;
	elseif ispl == 2
		isSPL = sessions.("63dB")==1;
	else
		isSPL = sessions.("83dB")==1;
	end

	passing_sesh = isBE & isSPL;
	pass_ind = find(passing_sesh);
	num_sessions = sum(passing_sesh);

	red = [187, 249, 186]/255;
	pink = [23, 64, 86]/255;
	color_gradient = [linspace(red(1),pink(1),num_sessions)', linspace(red(2),pink(2),num_sessions)',...
		linspace(red(3),pink(3),num_sessions)'];
	CFs = sessions.CF(pass_ind);
	[~, order] = sort(CFs);

	count = 0;
	ind = 1;

	for isesh = 1:num_sessions

		% Load in session
		session = sessions.Session{pass_ind(isesh)};
		TT = sessions.TT(pass_ind(isesh));
		N = sessions.N(pass_ind(isesh));
		CF = sessions.CF(pass_ind(isesh));
		filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
		load(fullfile(path, filename))

		[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
			= finddata(params);

		% Find dataset of interest
		if ispl == 1
			interest = find(bin & F0_200 & level_43);
		elseif ispl == 2
			interest = find(bin & F0_200 & level_63);
		elseif ispl == 3
			interest = find(bin & F0_200 & level_83);
		end
		if any(interest)
			data_ST = data{interest};
			rates = data_ST.rates_sm;
			fpeaks = data_ST.fpeaks;

			% Get CF
			%spont = population{isesh,1}{4, 1}.spont;

			% Get peak
			[~, peak_ind] = max(rates);
			peak_freq = fpeaks(peak_ind);

			% Align 
			fpeaks_re_CF = fpeaks-CF;
			fpeaks_re_max = fpeaks-peak_freq;

			% Normalize
			%rates = zscore(rates);
			rates = rates/max(rates);

			% Plot
			nexttile(ispl)
			hold on
			plot(fpeaks_re_max, rates, 'linewidth', 1.5, 'Color',color_gradient(order(isesh), :));
		end

	end
end

spls = [43, 63, 83];
for ispl = 1:3
	% Plotting Params
	nexttile(ispl)
	xlabel('Subtracted by CF (Hz)')
	ylabel('Avg. Rate (sp/s)')
	ylabel('z-score')
	box on
	title(['BE, ' num2str(spls(ispl)) ' dB SPL, F0=200, Bin'])
	set(gca, 'fontsize', 14)
end


%% Scatter plot of CF vs HHBW 

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BE');
num_sessions = sum(isBE);
BE_ind = find(isBE);
i = 1;
index = [];
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = (bin & F0_200 & level_43)|(bin & F0_200 & level_63)|(bin & F0_200 & level_83);
	if any(interest)
		index(i) = BE_ind(isesh);
		i = i+1;
	end
end

num_sessions = length(index);
width_CF = NaN(num_sessions, 3);
width = NaN(num_sessions, 3);
CFs = NaN(num_sessions);
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{index(isesh)};
	TT = sessions.TT(index(isesh));
	N = sessions.N(index(isesh));
	CF = sessions.CF(index(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest43 = find(bin & F0_200 & level_43);
	interest63 = find(bin & F0_200 & level_63);
	interest83 = find(bin & F0_200 & level_83);
	interest = [interest43', interest63',  interest83'];
	if any(interest)
		name = 'BE';
		CFs(isesh) = CF;
		num_interest = length(interest);
		for ind = 1:num_interest
			ids = interest(ind);
			data_ST = data{ids};
			param = params{ids};

			% Get CF rate
			[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
			CF_rate = data_ST.rates_sm(CF_ind);

			% Get spl 
			spl_shift = checkIfSPLShift(param.animal, param.session);
			spl = param.spl + spl_shift;	

			if spl == 43
				width(isesh, 1) = data_ST.width;
				if data_ST.width == 2400
					width_CF(isesh, 1) = 0;
				else
					width_CF(isesh, 1) = data_ST.width/CF;
				end
			elseif spl == 63
				width(isesh, 2) = data_ST.width;
				if data_ST.width == 2400
					width_CF(isesh, 2) = 0;
				else
					width_CF(isesh, 2) = data_ST.width/CF;
				end
			%elseif spl == 73
			%	width(isesh, 3) = data.width;
			elseif spl == 83
				width(isesh, 3) = data_ST.width;
				if data_ST.width == 2400
					width_CF(isesh, 3) = 0;
				else
					width_CF(isesh, 3) = data_ST.width/CF;
				end
			end			
		end
	end
end

% % Plot 
% width(width==0) = NaN;
% figure('position', [519,566,983,272])
% tiledlayout(1, 3, 'TileSpacing','tight', 'Padding','compact')
% spls = [43, 63, 83];
% for ispl = 1:3
% 
% 	nexttile
% 	scatter(CFs, width(:,ispl), 'filled')
% 	box on
% 	title([num2str(spls(ispl)) 'dB SPL'])
% 	xlabel('CF')
% 	ylabel('1/2 Height BW (Hz)')
% 	ylim([0 2400])
% 	set(gca, 'fontsize', 14)
% 	hold on
% 	set(gca, 'XScale', 'log')
% 	set(gca, 'YScale', 'log')
% 	xticks([0 200 500 1000 2000 5000 10000]/1000)
% 	yticks([0 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
% 
% 	% Add fit line
% 	single_width = width(:,ispl);
% 	nan_width = find(~isnan(single_width));
% 	single_width = single_width(nan_width);
% 	CF_single = CFs(nan_width);
% 	single_width = log10(single_width);
% 	CF_log = log10(CF_single);
% 
% 	p = polyfit(CF_log, single_width, 1); % Fit a polynomial of degree 1 (linear fit)
% 	x_fit = linspace(min(CF_single), max(CF_single), 100); % Generate x values for the fit line
% 	y_fit_log = polyval(p, x_fit); % Evaluate the polynomial at x_fit
% 	y_fit = exp(y_fit_log);
% 	plot(x_fit, y_fit_log, 'r-', 'DisplayName', 'Fit Line');
% 
% 	legend('Data', 'Fit Line', 'location', 'southeast')
% end

% Plot 
width(width==0) = NaN;
figure('position', [519,533,1150,305])
tiledlayout(1, 3, 'TileSpacing','tight', 'Padding','compact')
spls = [43, 63, 83];
for ispl = 1:3
	
	nexttile
	scatter(CFs/1000, width_CF(:,ispl), 'filled')
	box on
	number = width_CF(:,ispl);
	number(isnan(number)) = [];
	title([num2str(spls(ispl)) ' dB SPL, n=' num2str(length(number))])
	xlabel('CF (kHz)')
	if ispl == 1
		ylabel('1/2 Height BW / CF')
	end
	ylim([0.02 5])
	xlim([0.3 8])
	set(gca, 'fontsize', 17)
	hold on
	set(gca, 'XScale', 'log')
	set(gca, 'YScale', 'log')
	xticks([0 200 500 1000 2000 5000 10000]/1000)
	yticks([0 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
	yline(1, '--k', 'linewidth', 2)
	grid on

end

%% Plot imagesc of all BS responses sorted by CF

figure('position', [60,30,966,867])
tiledlayout(1, 3, 'Padding','compact', 'TileSpacing','none')
backgroundcolor =  [0.8 0.8 0.8];
spls = {'43', '63', '83'};

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BS');
num_sessions = sum(isBE);
BE_ind = find(isBE);

i = 1;
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = (bin & F0_200 & level_43)|(bin & F0_200 & level_63)|(bin & F0_200 & level_83);
	if any(interest)
		pass_ind(i) = BE_ind(isesh);
		i = i+1;
	end
end

%isSPL = sessions.("43dB")==1|sessions.("63dB")==1|sessions.("83dB")==1;
%isBE = strcmp(sessions.MTF, 'BS');
for ispl = 1:3

	% if ispl == 1
	% 	isSPL = sessions.("43dB")==1;
	% elseif ispl == 2
	% 	isSPL = sessions.("63dB")==1;
	% else
	% 	isSPL = sessions.("83dB")==1;
	% end

	%passing_sesh = isBE & isSPL;
	%pass_ind = find(passing_sesh);
	num_sessions = length(pass_ind);
	array_z = NaN(num_sessions,10000);
	CFs = NaN(num_sessions, 1);
	CF_names = cell(num_sessions, 1);
	for isesh = 1:num_sessions

		% Load in session
		session = sessions.Session{pass_ind(isesh)};
		TT = sessions.TT(pass_ind(isesh));
		N = sessions.N(pass_ind(isesh));
		CF = sessions.CF(pass_ind(isesh));
		filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
		load(fullfile(path, filename))

		[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
			= finddata(params);

		
		% Find dataset of interest
		if ispl == 1
			interest = find(bin & F0_200 & level_43);
		elseif ispl == 2
			interest = find(bin & F0_200 & level_63);
		else
			interest = find(bin & F0_200 & level_83);
		end

		CFs(isesh) = CF;
		if strcmp(session, 'R024S470')
			CF_names{isesh} = [num2str(round(CF)) ' Hz*'];
		elseif strcmp(session, 'R024S478')
			CF_names{isesh} = [num2str(round(CF)) ' Hz**'];
		else
			CF_names{isesh} = [num2str(round(CF)) ' Hz'];
		end
		num_interest = length(interest);
		for ind = 1:num_interest

			ids = interest(ind);
			data_ST = data{ids};
			param_ST = params{ids};
			spont = data{4}.spont;

			% General analysis
			rate = data_ST.rate;
			rate = rate - spont;
			fpeaks = data_ST.fpeaks;
			fpeaks_re_CF = log2(fpeaks/CF);
			%fpeaks_re_CF = fpeaks/CF;

			% Align by CF (approximately)
			f = linspace(-3, 3, 10000);
			%f = linspace(-1500, 1500, 100);
			[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
			[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
			f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

			% Interpolate & get z-score
			r_interp = interp1(fpeaks_re_CF, rate,f_interp, 'spline');
			z_rate = zscore(r_interp);
			%z_rate = r_interp;
			array_z(isesh, f_ind(1):f_ind(2)-1) = z_rate;
		end
	end

	% Order by CF 
	[~, max_ind] = sort(CFs);
	CF_order = array_z(max_ind,:);
	CFs_ordered = CF_names(max_ind);

	% Order by peak prominence 


	%nexttile
	s(ispl) = subplot(1, 3, ispl);
	h = imagesc(f, 1:size(CF_order, 1), CF_order);
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
	xlim([-2.5 2.5])
	clim([-4 4])
	xticklabels({'-2', '-1', 'CF', '1', '2'})
	for ind = 1:size(CF_order, 1)
		CF_label = CFs_ordered{ind};
		interval = 0.91 / size(CF_order, 1);
		start = 0.01;
		annotation('textbox',[start 0.935-interval*ind 0.07 0.049], ...
			'String',CF_label, 'EdgeColor',...
			'none', 'FontSize',12);
	end
	set(gca, 'fontsize', 14)
	title([spls{ispl} ' dB SPL'], 'fontsize', 18)
end

set(s(1), 'position', [0.07,0.0552, 0.305,0.9071])
set(s(2), 'position', [0.375,0.0552,0.305,0.9071])
set(s(3), 'position', [0.68,0.0552,0.305,0.9071])
colorbar

% BS 
annotation('textbox',[0.283 0.767 0.0866 0.0209],...
	'String','Neuron #1',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation('textbox',[0.586 0.767 0.0866 0.0209],...
	'String','Neuron #1',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation('textbox',[0.865 0.767 0.0866 0.0209],...
	'String','Neuron #1',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation('textbox',[0.2793 0.755 0.0866 0.0209],...
	'String','Neuron #2',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation('textbox',[0.579 0.755 0.0866 0.0209],...
	'String','Neuron #2',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation('textbox',[0.860 0.755 0.0866 0.0209],...
	'String','Neuron #2',...
	'FontSize',12,...
	'FitBoxToText','off',...
	'EdgeColor','none');

%% Plot imagesc of all BE responses sorted by CF

figure('position', [60,30,966,585])
tiledlayout(1, 3, 'Padding','compact', 'TileSpacing','none')
backgroundcolor =  [0.8 0.8 0.8];
spls = {'43', '63', '83'};

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BE');
num_sessions = sum(isBE);
BE_ind = find(isBE);
i = 1;
pass_ind = [];
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = (bin & F0_200 & level_43)|(bin & F0_200 & level_63)|(bin & F0_200 & level_83);
	if any(interest)
		pass_ind(i) = BE_ind(isesh);
		i = i+1;
	end
end

for ispl = 1:3

	num_sessions = length(pass_ind);
	array_z = NaN(num_sessions,10000);
	CFs = NaN(num_sessions, 1);
	CF_names = cell(num_sessions, 1);
	for isesh = 1:num_sessions

		% Load in session
		session = sessions.Session{pass_ind(isesh)};
		TT = sessions.TT(pass_ind(isesh));
		N = sessions.N(pass_ind(isesh));
		CF = sessions.CF(pass_ind(isesh));
		filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
		load(fullfile(path, filename))

		[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
			= finddata(params);

		% Find dataset of interest
		if ispl == 1
			interest = find(bin & F0_200 & level_43);
		elseif ispl == 2
			interest = find(bin & F0_200 & level_63);
		else
			interest = find(bin & F0_200 & level_83);
		end

		CFs(isesh) = CF;
		if strcmp(session, 'R027S247') & TT==3
			CF_names{isesh} = [num2str(round(CF)) ' Hz*'];
		elseif strcmp(session, 'R029S105')
			CF_names{isesh} = [num2str(round(CF)) ' Hz**'];
		else
			CF_names{isesh} = [num2str(round(CF)) ' Hz'];
		end
		num_interest = length(interest);
		for ind = 1:num_interest

			ids = interest(ind);
			data_ST = data{ids};
			param_ST = params{ids};
			spont = data{4}.spont;

			% General analysis
			rate = data_ST.rate;
			rate = rate - spont;
			fpeaks = data_ST.fpeaks;
			fpeaks_re_CF = log2(fpeaks/CF);
			%fpeaks_re_CF = fpeaks/CF;

			% Align by CF (approximately)
			f = linspace(-3, 3, 10000);
			%f = linspace(-1500, 1500, 100);
			[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
			[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
			f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

			% Interpolate & get z-score
			r_interp = interp1(fpeaks_re_CF, rate,f_interp, 'spline');
			z_rate = zscore(r_interp);
			%z_rate = r_interp;
			array_z(isesh, f_ind(1):f_ind(2)-1) = z_rate;
		end
	end

	% Order by CF 
	[~, max_ind] = sort(CFs);
	CF_order = array_z(max_ind,:);
	CFs_ordered = CF_names(max_ind);

	% Plot 
	s(ispl) = subplot(1, 3, ispl);
	h = imagesc(f, 1:size(CF_order, 1), CF_order);
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
	clim([-4 4])
	xticklabels({'-2', '-1', 'CF', '1', '2'})
	for ind = 1:size(CF_order, 1)
		CF_label = CFs_ordered{ind};
		interval = 0.873 / size(CF_order, 1);
		start = 0.01;
		annotation('textbox',[start 0.925-interval*ind 0.07 0.049], ...
			'String',CF_label, 'EdgeColor',...
			'none', 'FontSize',12);
	end
	set(gca, 'fontsize', 14)
	title([spls{ispl} ' dB SPL'], 'fontsize', 18)
end

set(s(1), 'position', [0.07,0.0752, 0.305,0.87])
set(s(2), 'position', [0.375,0.0752,0.305,0.87])
set(s(3), 'position', [0.68,0.0752,0.313,0.87])
colorbar

% BE Annotations
% Create textbox
annotation('textbox',...
	[0.600422832980975 0.881905982905982 0.0927589852008457 0.0452991452991452],...
	'String',{'Neuron #1'},...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.54545454545455 0.109256410256409 0.0927589852008457 0.0452991452991453],...
	'String','Neuron #2',...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.871504116712412 0.883615384615384 0.0927589852008457 0.0452991452991452],...
	'String',{'Neuron #1'},...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.824817402684954 0.110965811965811 0.0927589852008457 0.0452991452991453],...
	'String','Neuron #2',...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.29921999133324 0.881905982905982 0.0927589852008457 0.0452991452991452],...
	'String',{'Neuron #1'},...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

% Create textbox
annotation('textbox',...
	[0.241124227979641 0.110965811965811 0.0927589852008457 0.0452991452991453],...
	'String','Neuron #2',...
	'FontSize',14,...
	'FitBoxToText','off',...
	'EdgeColor','none');

%% Change in HHBW over level 

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BS');
num_sessions = sum(isBE);
BE_ind = find(isBE);
i = 1;
index = [];
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = (bin & F0_200 & level_43)|(bin & F0_200 & level_63)|(bin & F0_200 & level_83);
	if any(interest)
		index(i) = BE_ind(isesh);
		i = i+1;
	end
end

num_sessions = length(index);
width_bin = NaN(num_sessions, 3);
for iclus = 1:num_sessions

	% Load in data
	session = sessions.Session{index(iclus)};
	TT = sessions.TT(index(iclus));
	N = sessions.N(index(iclus));
	CF = sessions.CF(index(iclus));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even]...
		= finddata(params);

	for ispl = 1:3
		if ispl == 1
			ds_bin = find(bin & F0_200 & level_43);
			if any(ds_bin)
				width_bin(iclus, 1) = data{ds_bin(1)}.width/CF;
			end
		elseif ispl == 2
			ds_bin = find(bin & F0_200 & level_63);
			if any(ds_bin)
				width_bin(iclus, 2) = data{ds_bin(1)}.width/CF;
			end
		else
			ds_bin = find(bin & F0_200 & level_83);
			if any(ds_bin)
				width_bin(iclus, 3) = data{ds_bin(1)}.width/CF;
			end
		end
	end
end
%%
figure('position', [890,461,775,305])
tiledlayout(1, 2, 'Padding','compact', 'TileSpacing','compact')

% nexttile
% % edges = linspace(-1500, 1500, 21);
% edges = linspace(-1, 1, 21);
% histogram(differences, edges)
% hold on
% xline(0, 'k', 'LineWidth',1.5)
% xline(mean(differences, 'omitnan'), 'r', 'LineWidth',1.5)
% xline(median(differences, 'omitnan'), '--r', 'LineWidth',1.5)
% set(gca, 'fontsize', 18)
% title('1/2 Height Differences')
% xlabel('83dB SPL - 43 dB SPL (Hz)')
% ylabel('# Neurons')

differences = NaN(length(width_bin), 1);
for ii = 1:length(width_bin)
	if ~isnan(width_bin(ii,3)) && ~isnan(width_bin(ii,1))
		differences(ii) = width_bin((ii),3)-width_bin((ii),1);
	elseif ~isnan(width_bin(ii,2)) && ~isnan(width_bin(ii,1))
		differences(ii) = width_bin((ii),2)-width_bin((ii),1);
	elseif ~isnan(width_bin(ii,3)) && ~isnan(width_bin(ii,2))
		differences(ii) = width_bin((ii),3)-width_bin((ii),2);
	end
end
decrease = find(sign(differences)==-1);
increase = find(sign(differences)==1);
% steady = abs(differences)<=0.1;
% decrease = differences<0.1;
% increase = differences>0.1;

% nexttile
% hold on
% plot([43, 63, 83], width_bin(steady,:)', 'color', 'k', 'LineWidth',1.5)
% xticks([43, 63, 83])
% set(gca, 'fontsize', 18)
% title(['Steady Peak (n=' num2str(sum(steady)) ')'])
% xlabel('Level (dB SPL)')
% ylabel('1/2 Height BW / CF (Hz)')
% box on
% xlim([40 86])


nexttile
hold on
plot([43, 63, 83], width_bin(increase,:)', 'color', [0.85,0.11,0.38, 0.6], 'LineWidth',1.5)
xticks([43, 63, 83])
set(gca, 'fontsize', 18)
title(['Broadening Peak (n=' num2str(length(increase)) ')'])
xlabel('Level (dB SPL)')
ylabel('1/2 Height BW / CF (Hz)')
box on
xlim([40 86])
plot([43, 63, 83], mean(width_bin(increase,:), 'omitnan'), 'k', 'LineWidth',2)
plot([43, 63, 83], median(width_bin(increase,:), 'omitnan'), '--k', 'LineWidth',2)
ylim([0 2])

nexttile
hold on
plot([43, 63, 83], width_bin(decrease,:)', 'color',[0.12,0.53,0.90, 0.6], 'LineWidth',1.5)
xticks([43, 63, 83])
set(gca, 'fontsize', 18)
title(['Sharpening Peak (n=' num2str(length(decrease)) ')'])
xlabel('Level (dB SPL)')
box on
xlim([40 86])
ylim([0 2])
yticklabels([])
plot([43, 63, 83], mean(width_bin(decrease,:), 'omitnan'), 'k', 'LineWidth',2)
plot([43, 63, 83], median(width_bin(decrease,:), 'omitnan'), '--k', 'LineWidth',2)
legend('', '', '', '', '', '', '', '', '', '','','','','','','', '', '','Mean', 'Median', 'fontsize', 18)


%% Energy matched to BE responses

% Find BE sessions 
isBE = strcmp(sessions.MTF, 'BE');
num_sessions = sum(isBE);
BE_ind = find(isBE);
i = 1;
index = [];
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{BE_ind(isesh)};
	TT = sessions.TT(BE_ind(isesh));
	N = sessions.N(BE_ind(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest = (bin & F0_200 & level_43)|(bin & F0_200 & level_63)|(bin & F0_200 & level_83);
	if any(interest)
		index(i) = BE_ind(isesh);
		i = i+1;
	end
end

num_sessions = length(index);
r2_all = NaN(num_sessions, 3);
CFs = NaN(num_sessions);
for isesh = 1:num_sessions

	% Load in session
	session = sessions.Session{index(isesh)};
	TT = sessions.TT(index(isesh));
	N = sessions.N(index(isesh));
	CF = sessions.CF(index(isesh));
	filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
	load(fullfile(path, filename))

	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83]...
		= finddata(params);

	% Find dataset of interest
	interest43 = find(bin & F0_200 & level_43);
	interest63 = find(bin & F0_200 & level_63);
	interest83 = find(bin & F0_200 & level_83);
	interest = [interest43', interest63',  interest83'];
	if any(interest)
		name = 'BE';
		CFs(isesh) = CF;
		num_interest = length(interest);
		for ind = 1:num_interest
			ids = interest(ind);
			data_ST = data{ids};
			param = params{ids};

			% Get CF rate
			[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
			CF_rate = data_ST.rates_sm(CF_ind);

			% Get spl 
			spl_shift = checkIfSPLShift(param.animal, param.session);
			spl = param.spl + spl_shift;

			% Create stimulus
			param.Fs = 100000;
			param.physio = 1;
			param.mnrep = 20;
			param.dur = 0.3;
			[param] = generate_synthetictimbre(param);
			param.num_stim = size(param.stim, 1);

			% Energy model
			Fs = 100000;
			stimulus = [param.stim zeros(size(param.stim,1),0.1*Fs)];
			gamma_param.srate = Fs;
			tvals = (1:length(stimulus))/Fs;
			gamma_IF_reg = zeros(length(CFs),length(tvals));
			impaired = 0; % 0 = not impaired; 1 = 'impaired'

			pin_gamma = zeros(size(stimulus, 1), Fs*param.dur+0.1*Fs);
			for istim = 1:size(stimulus, 1)
				gamma_param.fc = CF;
				pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
			end
			pin_gamma = pin_gamma(:,1:param.dur*Fs);
			energy = sqrt(mean(pin_gamma.^2,2));

			% R^2
			[rate, ~] = plotSyntheticTimbre(param, energy, 0);
			R_int = corrcoef(data_ST.rate,rate);
			r2 = R_int(1, 2).^2;

			if spl == 43
				r2_all(isesh, 1) = r2;
			elseif spl == 63
				r2_all(isesh, 2) = r2;
			elseif spl == 83
				r2_all(isesh, 3) = r2;
			end			
		end
	end
end


%% Plot 
figure('position', [519,592,1100,240])
tiledlayout(1, 3, 'TileSpacing','tight', 'Padding','compact')
spls = [43, 63, 83];
for ispl = 1:3
	
	nexttile
	edges = linspace(0, 1, 21);
	histogram(r2_all(:,ispl), edges)
	box on
	number = r2_all(:,ispl);
	number(isnan(number)) = [];
	title([num2str(spls(ispl)) ' dB SPL, n=' num2str(length(number))])
	xlabel('Energy Prediction R^2')
	if ispl == 1
		ylabel('# Neurons')
	end
	%ylim([0.02 5])
	xlim([0 1])
	set(gca, 'fontsize', 16)
	hold on
	grid on

end


%% Functions 

function [bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even] = finddata(population)

	% Find binaural
	bin = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==2, population);

	% Find contra
	contra = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==1, population);

	% Find F0 = 100Hz
	F0_100 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==100, population);

	% Find F0 = 200Hz
	F0_200 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==200, population);

	% Find 43 db SPL
	level_43 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==43 || p.spl==40), population);

	% Find 63 dB SPL
	level_63 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==63 || p.spl==60), population);

	% Find 73 dB SPL
	level_73 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==73 || p.spl==70), population);

	% Find 83 dB SPL
	level_83 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==83 || p.spl==80), population);

	% Find even 200Hz fpeak_mid
	even = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);

end