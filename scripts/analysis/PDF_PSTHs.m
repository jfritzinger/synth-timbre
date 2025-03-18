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
filename = 'PSTH_VS_all';
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

% Find sessions with MTF
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);
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
		tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','tight', ...
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
		title(['MTF: ' MTF_target], 'FontSize', fontsize)
		ylabel('Avg. Rate (sp/s)')
		%legend('Unmodulated', 'Location','southwest')

		% Plot synthetic timbre
		ileg = 1;
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		ispl = 2;
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
			leg{ileg} = '';
			ileg = ileg + 1;

			% Plot BE Synth Timbre
			nexttile(3)
			hold on
			plot(param.fpeaks, temporal.VS)
			num_fpeaks = length(param.fpeaks);
			plot(param.fpeaks, smooth_rates(temporal.VS,zeros(num_fpeaks, 1),...
				ones(num_fpeaks, 1), CF), 'k')
		end
		plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
		xline(CF, '--', 'Color',CF_color, 'linewidth', 1.5)
		xlim(plot_range)
		leg{ileg} = 'CF';
		xline(CF, '--')
		xlabel('Spectral Peak Freq. (Hz)')
		ylabel('Vector Strength')
		title('Example neuron VS, BS')
		grid on
		ylim([0 1])
		set(gca, 'fontsize', fontsize)
		legend(leg, 'Location','best')
		clear leg

		% Plot PSTH
		freq_values = round(param.fpeaks);
		if bin200_MTF(ineuron, ispl)==1
			nexttile([3, 1])
			hold on
			max_rate = max(temporal.PSTH, [], 'all');
			for j = 1:num_fpeaks

				% Plot PSTHs
				counts = temporal.PSTH(j,:);
				edges = temporal.t;
				t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
				x_patch = repelem(edges, 2);
				y_patch = repelem([0; counts(:); 0]', 2);
				y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
				offset = (j-1)*max_rate; % Adjust offset amount
				patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');

				% Plot smoothed PSTH
				plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);

			end
			ylim([0 max_rate*num_fpeaks])
			xlabel('Time (ms)')
			box on
			yticks(linspace(max_rate/2, max_rate*num_fpeaks-max_rate/2, num_fpeaks))
			yticklabels(freq_values)
			title('PSTH')
			set(gca, 'fontsize', fontsize)
		end

		% Plot period histogram
		if bin200_MTF(ineuron, ispl)==1
			nexttile([3, 1])
			max_rate = max(temporal.p_hist, [], 'all');
			for j = 1:num_fpeaks

				% Plot PSTHs
				counts = temporal.p_hist(j,:);
				edges = temporal.t_hist;
				t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
				x_patch = repelem(edges, 2);
				y_patch = repelem([0; counts(:); 0]', 2);
				y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
				offset = (j-1)*max_rate; % Adjust offset amount
				patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
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
		end

		% Add to PDF
		[plt1, images] = addtoSTPDF(images, fig, putative);
		append(rpt, plt1);

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
