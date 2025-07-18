%% PDF of all synthetic timbre responses (binaural or contra)

clear
import mlreportgen.dom.*
import mlreportgen.report.*

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', ...
	spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
filename = 'Period_Hist_All';
images = {}; %hold all plots as images, need to delete when finished
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss')
report_name = sprintf('%s/pdfs/%s_%s.pdf', savepath, datetime, filename);
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.01in';
pm.PageMargins.Header = '0.01in';
pm.PageMargins.Bottom = '0.01in';
pm.PageMargins.Footer = '0.01in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

%% Plot each dataset
warning('off')

% Find sessions with MTF
%MTF_target = 'BS';
%isMTF = strcmp(sessions.MTF, MTF_target);
%isMTF = contains(sessions.MTF, MTF_target);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
bin200_MTF = bin200; % & isMTF;
has_data = bin200_MTF(:,1) | bin200_MTF(:,2) | bin200_MTF(:,3) | bin200_MTF(:,4);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Plot each neuron
spls = [43, 63, 73, 83];
for isesh = 1:num_sessions
	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))

		% Paragraph intro
		label = sprintf("%s, CF = %0.0fHz, %s\n", putative, CF, MTF_shape);
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		% Output
		fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);

		% Set up figure
		fig = figure('Position',[292,274,1264,420]);
		tiledlayout(4, 3, 'TileSpacing','compact', 'Padding','tight', ...
			'TileIndexing','columnmajor')
		x_label = [1000 2000 3000 4000 6000 8000]./1000;
		fontsize = 10;

		% Plot RM
		params_RM = data{2, 2};
		if ~isempty(params_RM)
			data_RM = analyzeRM(params_RM);
			nexttile(1)
			hold on
			spont_color = [0.4 0.4 0.4];
			CF_color = [0.3 0.3 0.3];
			plot(data_RM.freqs./1000,data_RM.rates(:,5),'color', '#20116B','LineWidth',2) % 73 dB
			plot(data_RM.freqs./1000,data_RM.rates(:,4),'color', '#5E50A9','LineWidth',2) % 53 dB
			plot(data_RM.freqs./1000,data_RM.rates(:,3),'color', '#A49BD0','LineWidth',2) % 33 dB
			plot(data_RM.freqs([1 end])./1000,[1 1]*data_RM.spont,'-','LineWidth',2, 'Color',spont_color)
			xline(CF/1000, '--', 'Color', CF_color,'LineWidth',2);
			box on
			grid on
			hold off
			ylim([0 max(data_RM.rates, [], 'all')+10])
			set(gca,'XTick',[])
			xlim([250 14000]./1000)
			xticks(x_label)
			set(gca, 'XScale', 'log');
			set(gcf, 'color', 'w')
			set(gca,'fontsize',fontsize)
			ylabel('Avg. Rate (sp/s)')
			%legend('73dB SPL', '53dB SPL', '33dB SPL', 'Spont. Rate', 'fontsize',fontsize-2, 'location', 'best')
			title('Response Map')
			xlabel('Frequency (kHz)')
		end

		% Plot MTF
		params_MTF = data{3, 2};
		data_MTF = analyzeMTF(params_MTF);
		nexttile(2)
		hold on
		line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',2);
		errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',2, 'Color','k');
		line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 2);
		hold off
		set(gca, 'XScale', 'log');
		xlim([data_MTF.fms(1) data_MTF.fms(end)])
		xticks([1 2 5 10 20 50 100 200 500])
		%xlabel('Modulation Frequency')
		%ylim(MTF_range)
		set(gca,'fontsize',fontsize)
		ylimit = ylim;
		ylim([0 ylimit(2)])
		grid on
		box on
		%title(['MTF: ' MTF_target], 'FontSize', fontsize)
		ylabel('Avg. Rate (sp/s)')
		%legend('Unmodulated', 'Location','southwest')

		% Plot synth timbre average rates
		ileg = 1;
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		nexttile
		for ispl = 1:4
			if bin200_MTF(ineuron, ispl)==1
				param_ST = data(5+ispl, 2);
				data_ST = analyzeST(param_ST, CF);
				data_ST = data_ST{1};

				% Get CF rate
				[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
				CF_rate = data_ST.rates_sm(CF_ind);

				% Get spl
				leg{ileg} = [num2str(param_ST{1}.spl) ' dB SPL'];
				ileg = ileg + 1;

				% Plot BE Synth Timbre
				
				hold on
				plot(data_ST.fpeaks,data_ST.rates_sm, 'linewidth', 1.5,...
					'Color',data_colors{ispl});
			end
		end
		plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
		xline(CF, '--', 'Color',CF_color, 'linewidth', 1.5)
		xlim(plot_range)	
		xlabel('Frequency (Hz)')
		ylabel('Avg. Rate (sp/s)')
		box on
		title('Synth Timbre: F0=200, Bin')
		grid on
		set(gca, 'fontsize', fontsize)

		% Plot synthetic timbre
		ileg = 1;
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		nexttile
		clear leg
		for ispl = 1:4
			if bin200_MTF(ineuron, ispl)==1
				param_ST = data(5+ispl, 2);
				data_ST = analyzeST(param_ST, CF);
				data_ST = data_ST{1};
				param = param_ST{1};
				temporal = analyzeST_Temporal(param, data_ST);

				% Get CF rate
				[~, CF_ind] = min(abs(CF-data_ST.fpeaks));
				CF_rate = data_ST.rates_sm(CF_ind);

				% Get spl
				leg{ileg} = [num2str(param_ST{1}.spl) ' dB SPL'];
				ileg = ileg + 1;
				% leg{ileg} = '';
				% ileg = ileg + 1;

				% Plot BE Synth Timbre
				hold on
				plot(param.fpeaks, temporal.VS, 'Color',data_colors{ispl})
				% num_fpeaks = length(param.fpeaks);
				% plot(param.fpeaks, smooth_rates(temporal.VS,zeros(num_fpeaks, 1),...
				% 	ones(num_fpeaks, 1), CF), 'k')
			end
		end
		plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
		xline(CF, '--', 'Color',CF_color, 'linewidth', 1.5)
		xlim(plot_range)
		leg{ileg} = 'CF';
		xline(CF, '--')
		xlabel('Spectral Peak Freq. (Hz)')
		ylabel('Vector Strength')
		title('Vector Strength')
		grid on
		ylim([0 1])
		set(gca, 'fontsize', fontsize)
		legend(leg, 'Location','best')
		clear leg

		% Plot PSTH
		freq_values = round(param.fpeaks);
		% if bin200_MTF(ineuron, ispl)==1
		% 	nexttile([3, 1])
		% 	hold on
		% 	max_rate = max(temporal.PSTH, [], 'all');
		% 	for j = 1:num_fpeaks
		% 
			% 	% Plot PSTHs
			% 	counts = temporal.PSTH(j,:);
			% 	edges = temporal.t;
			% 	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
			% 	x_patch = repelem(edges, 2);
			% 	y_patch = repelem([0; counts(:); 0]', 2);
			% 	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
			% 	offset = (j-1)*max_rate; % Adjust offset amount
			% 	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
		% 
			% 	% Plot smoothed PSTH
			% 	plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);
		% 
		% 	end
		% 	ylim([0 max_rate*num_fpeaks])
		% 	xlabel('Time (ms)')
		% 	box on
		% 	yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
		% 	yticklabels(freq_values)
		% 	title('PSTH')
		% 	set(gca, 'fontsize', fontsize)
		% end

		% Plot period histogram
		hist_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;...
			0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
		nexttile([4, 1])
		hold on
		max_rate = max(temporal.p_hist, [], 'all');
		for ispl = 2
			if bin200_MTF(ineuron, ispl)==1
				param_ST = data(5+ispl, 2);
				data_ST = analyzeST(param_ST, CF);
				data_ST = data_ST{1};
				param = param_ST{1};
				num_fpeaks = length(param.fpeaks);
				temporal = analyzeST_Temporal(param, data_ST);			
				for j = 1:num_fpeaks
					counts = temporal.p_hist(j,:);
					edges = temporal.t_hist;
					t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
					x_patch = repelem(edges, 2);
					y_patch = repelem([0; counts(:); 0]', 2);
					y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
					offset = (j-1)*max_rate; % Adjust offset amount
					patch(x_patch, y_patch + offset, hist_colors(ispl, :),...
						'FaceAlpha',0.4, 'EdgeColor','k');
				end
			end
		end
		ylim([0 max_rate*num_fpeaks])
		xlabel('Time (ms)')
		box on
		yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
		yticklabels(freq_values)
		xlim([0 5])
		xticks(1:5)
		grid on
		title('Period Histogram')
		set(gca, 'fontsize', fontsize)

		% Plot phase locking to harmonics 
		nexttile([4 1])
		if bin200_MTF(ineuron, ispl)==1
			VS_harms2 =flipud(temporal.VS_harms);
			p_vals =flipud(temporal.VS_p_harms);
			sig = p_vals<0.01;
			VS_harms2(~sig) = NaN;
			peak_harm = param.fpeak_mid;
			imagesc(1:10, 1:40, VS_harms2)
			hold on
			for j = 1:num_fpeaks
				rectangle('position', [peak_harm-0.5 j-0.5, 1, 1], ...
					'EdgeColor','k', 'LineWidth',1.5)
			end
			xlim([0.51 10.51])
			xlabel('Harmonic Number')
			yticklabels([])
			clim([0 1])
			colorbar
		end

		% Add to PDF
		[plt1, images] = addtoSTPDF(images, fig, putative);
		append(rpt, plt1);

		% % Plot heatmap for each
		% for ispl = 1:4	
		% 	nexttile
		% 	if bin200_MTF(ineuron, ispl)==1
		% 
		% 		param_ST = data(5+ispl, 2);
		% 		data_ST = analyzeST(param_ST, CF);
		% 		data_ST = data_ST{1};
		% 		param = param_ST{1};
		% 
		% 		temporal = analyzeST_Temporal(param, data_ST);
		% 		num_fpeaks = length(data_ST.fpeaks);
		% 		max_rate = max(temporal.p_hist, [], 'all');
		% 		freq_values = round(param.fpeaks);
		% 		p_hist = temporal.p_hist;
		% 		edges = temporal.t_hist;
		% 		t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
		% 
		% 		grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
		% 		grayMap = flipud(grayMap);
		% 		h = pcolor(t_bin, data_ST.fpeaks, p_hist);
		% 		set(h, 'EdgeColor', 'none');
		% 		hold on
		% 		yline(CF, 'r', 'LineWidth',1)
		% 		title(sprintf('%d dB SPL, Period Hist', spls(ispl)));
		% 		colorbar;
		% 		%shading interp
		% 		axis square;
		% 		colormap(grayMap);
		% 		%xlim([0 5])
		% 		%clim([0 70])
		% 		ylabel('Spectral Peak Freq. (Hz)')
		% 		xlabel('Period (ms)')
		% 	end
		% end

		% % Calculate ISI for each rep
		% nreps = param.nrep;
		% ISI_all = cell(1, num_fpeaks);
		% nbins = 30;
		% edges = linspace(0, 20, nbins+1); % Bin edges
		% counts_all = zeros(num_fpeaks, nbins); % Pre-allocate counts_all
		% for jj = 1:num_fpeaks
		% 	x = temporal.x{jj} / 1000; % ms
		% 	y = temporal.y{jj};
		% 
		% 	ISI = arrayfun(@(ii) diff(x(y == ii)), 1:nreps, 'UniformOutput', false);
		% 	ISI_all{jj} = vertcat(ISI{:});
		% 
		% 	counts_all(jj, :) = histcounts(ISI_all{jj}, edges);
		% end
		% 
		% % Plot ISI histograms
		% nexttile([4,1])
		% max_rate = max(counts_all, [], 'all');
		% freq_values = round(param.fpeaks);
		% for j = 1:num_fpeaks
		% 	counts = counts_all(j,:);
		% 
		% 	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
		% 	x_patch = repelem(edges, 2);
		% 	y_patch = repelem([0; counts(:); 0]', 2);
		% 	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
		% 	offset = (j-1)*max_rate; % Adjust offset amount
		% 	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
		% end
		% ylim([0 max_rate*num_fpeaks])
		% xlabel('Time (ms)')
		% box on
		% yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
		% yticklabels(freq_values)
		% grid on
		% title('ISI Histogram')
		% set(gca, 'fontsize', fontsize)

		% %% Plots FFT
		% 
		% temporal = analyzeST_Temporal(param, data_ST);
		% num_fpeaks = length(data_ST.fpeaks);
		% spike_hist = temporal.PSTH;
		% edges = temporal.t;
		% bin_width = (edges(2) - edges(1))/1000; % seconds
		% fs = round(1/bin_width);
		% [VS, fft_output, f] = calcFFT(spike_hist, param, fs);
		% nexttile([4, 1])
		% max_rate = max(fft_output, [], 'all');
		% freq_values = round(param.fpeaks);
		% hold on
		% for j = 1:num_fpeaks % Plot FFTs
		% 	offset = (j-1)*max_rate; % Adjust offset amount
		% 	plot(f, fft_output(j,:) + offset);
		% end
		% ylim([0 max_rate*num_fpeaks])
		% xlabel('Frequency (Hz)')
		% box on
		% yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
		% yticklabels(freq_values)
		% grid on
		% title('FFT(PSTH)')
		% xlim([0 fs/2])
		% set(gca, 'fontsize', fontsize)


		% %% Plot ISI scatter plots
		% % Calculate ISI for each rep
		% ISI_all = cell(1, num_fpeaks);
		% counts_all = zeros(num_fpeaks, 50); % Pre-allocate counts_all
		% for jj = 1:num_fpeaks
		% 	ISI = [];
		% 	x_isi = [];
		% 	y_isi = [];
		% 	x = temporal.x{jj}/1000; % ms
		% 	y = temporal.y{jj};
		% 	for ii = 1:30
			% 	valid = y==ii;
			% 	isi = diff(x(valid));
			% 	ISI = [ISI isi'];
			% 	x_isi = [x_isi isi(1:end-1)'];
			% 	y_isi = [y_isi isi(2:end)'];
		% 	end
		% 	ISI_all{jj} = ISI;
		% 	x_isi_all{jj} = x_isi;
		% 	y_isi_all{jj} = y_isi;
		% end
		% 
		% 
		% fig = figure('Position',[560,42,1050,806]);
		% if num_fpeaks < 43
		% 	tiledlayout(6, 7, 'Padding','compact', 'TileSpacing','tight')
		% 	last_row = 36:41;
		% elseif num_fpeaks < 50
		% 	tiledlayout(7, 7, 'Padding','compact', 'TileSpacing','tight')
		% 	last_row = 42:49;
		% else
		% 	tiledlayout(8, 7, 'Padding','compact', 'TileSpacing','tight')
		% 	last_row = 50:56;
		% end
		% for jj = 1:num_fpeaks
		% 	T = 5; % period
		% 	nexttile
		% 	hold on
		% 	line(repmat([0;10*T],1,10),[1;1]*(1:10)*T,'Color','r', 'linewidth', 0.3)
		% 	line([1;1]*(1:10)*T,repmat([0;10*T],1,10),'Color','r', 'linewidth', 0.3)
		% 	axis square
		% 	scatter(x_isi_all{jj}, y_isi_all{jj},8, 'filled', 'k')
		% 
		% 	if ismember(jj, last_row)
			% 	xlabel('ISI n (ms)')
		% 	else
			% 	xticklabels([])
		% 	end
		% 
		% 	if ismember(jj, 1:7:41)
			% 	ylabel('ISI n+1 (ms)')
		% 	else
			% 	yticklabels([])
		% 	end
		% 	title(sprintf('%d Hz', freq_values(jj)))
		% 	xlim([0 30])
		% 	ylim([0 30])
		% end
		% 
		% % Add to PDF
		% [plt1, images] = addtoSTPDF(images, fig, [putative '_ISI']);
		% append(rpt, plt1);

	end
end

% Closes and opens PDF to view
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)

%% FUNCTIONS

function [img, images] = addtoSTPDF(images, fig, title)
import mlreportgen.dom.*

% Set figure size, recommended
values = [8, 10];
fig.PaperSize = values;
fig.PaperPosition = [0 0 values];
fig.Units = 'inches';
fig.Position(3:4) = values;

% Add the plot to the document
name = sprintf('%s.svg', title);
print(fig, name, '-dsvg');
img = Image(name);
delete(fig) %delete plot figure window
images = [images {img}];

end

%% FUNCTIONS 

function [VS, fft_output, f] = calcFFT(spike_hist, params, fs)

freq = 200;
for ii = 1:length(params.fpeaks)

	% Cut to integer # of cycles (50 cycles, excluding onset)
	onsetwin = 0.05; % ms
	spike_wo_onset = spike_hist(ii, round(onsetwin*fs):end);

	% Find frequency index for 200 Hz
	L = length(spike_wo_onset);

	% Take FFT 
	Y = fft(spike_wo_onset);
	P2 = abs(Y)/L; % Normalize by signal length
	P1 = P2(1:L/2+1);
	P1(2:end-1) = 2*P1(2:end-1);
	f = fs/L*(0:(L/2));

	% Get FFT at 200 Hz 
	ind = f==freq;
	%VS(ii) = P1(ind);
	fft_output(ii,:) = P1;
	VS = 0;

end
end