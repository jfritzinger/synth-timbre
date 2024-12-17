%% plot_binmode_example


%% Load in spreadsheet and set up PDF
clear

import mlreportgen.dom.*
import mlreportgen.report.*

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
report_path = 'figures/pdfs/';
filename = 'F0Comparison_All';
images = {}; %hold all plots as images, need to delete when finished
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss')
report_name = sprintf('%s%s_%s.pdf', fullfile(base, report_path), datetime, filename);
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

%% Determine sessions of interest

% Find sessions for target synthetic timbre response
% bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
% bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
% bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
% bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
% contra200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
% contra200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
% contra200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
% contra200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);
% both200 = bin200 & contra200;
% has_data = both200(:,1) | both200(:,2) | both200(:,3) | both200(:,4);
% index = find(has_data);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
contra200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_100);
contra200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_100);
contra200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_100);
both200 = bin200 & contra200;
has_data = both200(:,1) | both200(:,2) | both200(:,3);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

%% Plot
bin_colors = {'#3F985C', '#3F985C', '#3F985C', '#3F985C'};
contra_colors = {'#1267B9', '#1267B9', '#1267B9', '#1267B9'};
linewidth = 1.5;

for isesh = 1:num_sessions

	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, [putative '.mat']))

		% Paragraph intro
		% label = '';
		% p = Paragraph(label);
		% p.FontSize = "14pt";
		% p.WhiteSpace = "preserve";
		% append(rpt,p);

		label = sprintf("%s, CF = %0.0fHz, %s\n", putative, CF, MTF_shape);
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		% Output
		fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);


		fig = figure('Position',[70,607,951,250]);
		tiledlayout(1, 4, 'TileSpacing','compact', 'Padding','tight')
		for ind = 1:3 %num_ds

			% Plot
			nexttile
			hold on
			if both200(ineuron, ind)==1
				for ibin = 1:2
					if ibin == 1
						colors = contra_colors;
					else
						colors = bin_colors;
					end

					% Analysis (for bin and contra)
					% params = data(5+ind, 1:2);
					% data_ST  = analyzeST(params);
					% rate = data_ST{ibin}.rate;
					% rate_std = data_ST{ibin}.rate_std;
					% rlb = data_ST{ibin}.rlb;
					% rub = data_ST{ibin}.rub;
					% fpeaks = data_ST{ibin}.fpeaks;
					% spl = data_ST{ibin}.spl;
					% rate_sm = data_ST{ibin}.rates_sm;

					if ibin == 1 % 200 Hz
						if ind == 1 || ind == 2
							params = data(5+ind, 2);
						else
							params = data(5+ind+1, 2);
						end
					else % 100 Hz
						params = data(9+ind, 2);
					end
					data_ST  = analyzeST(params);
					rate = data_ST{1}.rate;
					rate_std = data_ST{1}.rate_std;
					rlb = data_ST{1}.rlb;
					rub = data_ST{1}.rub;
					fpeaks = data_ST{1}.fpeaks;
					spl = data_ST{1}.spl;
					rate_sm = data_ST{1}.rates_sm;

					% Plot
					rates_sm = smooth_rates(rate, rlb, rub);
					errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
						'linestyle', 'none', 'linewidth', 0.8, 'color', colors{ind})
					plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',colors{ind})
					%label{1} = [num2str(params{order(ind)}.spl) ' dB SPL'];

					% [~, peak_ind] = max(data_ST{ibin}.rate);
					% peak_f(ibin) = data_ST{ibin}.fpeaks(peak_ind);
					[~, peak_ind] = max(data_ST{1}.rate);
					peak_f(ibin) = data_ST{1}.fpeaks(peak_ind);
					width(ibin) = data_ST{1}.width;
				end

				plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)]./1000;
				xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
				set(gca, 'Fontsize', 6, 'XTick', plot_range(1)+0.200:0.400:plot_range(2)-0.200);
				xlim(plot_range);
				grid on
				ylimit = ylim;
				ylim([0 ylimit(2)])
				ylabel('Avg. rate (sp/s)')
				xlabel('Spectral Peak  Frequency (Hz)')
				title(sprintf('%d dB SPL', data_ST{1}.spl))
				box on
				if ind == 1
					legend({'200 Hz', '', '100 Hz'}, 'Location','best')
				end
				q_factor(ind,1:2) = [width(1)/peak_f(1) width(2)/peak_f(2)];
				q_factor2(ind,1:2) = [peak_f(1)/width(1) peak_f(2)/width(2)];

			end
		end

		% Q-value
		nexttile
		plot(q_factor2', 'o-')
		title('Q_H_H_B_W')
		ylabel('Freq_p_e_a_k / HHBW')
		xticks([1, 2])
		xlim([0.5 2.5])
		%xticklabels({'Contra', 'Bin'})
		xticklabels({'200 Hz', '100 Hz'})
		xlabel('SPL')
		set(gca, 'Fontsize', 6)
		legend('43', '63', '83', 'Location','best')
		clear q_factor2
		grid on

		% Add to PDF
		[plt1, images] = addtoSTPDF(images, fig, putative);
		append(rpt, plt1);

		label = '';
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

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
values = [8, 2.2];
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
