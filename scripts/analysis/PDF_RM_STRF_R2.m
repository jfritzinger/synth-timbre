%% plot_STRF_and_RM_R2
%
%
%
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
filename = 'RM_STRF_Predictions';
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


%% Analyze and plot data 

binmodes = {'Contra', 'Bin'};
SPLs = {'63', '73'};
for ispl = 1:2
	for ibin = 2

		if ispl == 1 && ibin == 1
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con) & cellfun(@(s) contains(s, 'R'), sessions.STRF_con);
		elseif ispl == 1 && ibin == 2
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
		elseif ispl == 2 && ibin == 1
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con) & cellfun(@(s) contains(s, 'R'), sessions.STRF_con);
		elseif ispl == 2 && ibin == 2
			has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
		end
		index = find(has_data);

		for isesh = 1:length(index)

			% Load in data
			s_ind = index(isesh);
			putative = sessions.Putative_Units{s_ind};
			CF = sessions.CF(s_ind);
			load(fullfile(datapath, 'neural_data', [putative '.mat']), 'data');
			params_STRF = data{4,ibin};
			param_ST = data(6+ispl,ibin); 
			param_RM = data{2, 2};

			% Load in STRF model
			if ispl == 1 && ibin == 1
				filename = [putative '_STRF_63_Contra'];
				load(fullfile(datapath, 'STRF_Models', [filename '.mat']), "STRFmodel_con")
			elseif ispl == 1 && ibin == 2
				filename = [putative '_STRF_63_Bin'];
				load(fullfile(datapath, 'STRF_Models', [filename '.mat']), "STRFmodel")
			elseif ispl == 2 && ibin == 1
				filename = [putative '_STRF_73_Contra'];
				load(fullfile(datapath, 'STRF_Models', [filename '.mat']), "STRFmodel_con")
			elseif ispl == 2 && ibin == 2
				filename = [putative '_STRF_73_Bin'];
				load(fullfile(datapath, 'STRF_Models', [filename '.mat']), "STRFmodel")
			end

			% Load in RM analysis 
			filename = sprintf('%s_RM_%s_%s', putative, SPLs{ispl}, binmodes{ibin});
			load(fullfile(datapath, 'RM_predictions', [filename '.mat']), "RM_R2")

			% General analysis
			data_STRF = analyzeSTRF(params_STRF);
			data_ST = analyzeST(param_ST);
			param_ST = param_ST{1};
			data_ST = data_ST{1};

			%% Plot 
			fig = figure('Position',[3,632,623,273]);
			tiledlayout(1, 3)

			params_MTF = data{3, 2};
			if ~isempty(params_MTF)
				data_MTF = analyzeMTF(params_MTF);
				nexttile
				hold on
				spont_color = [0.4 0.4 0.4];
				line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',1.5);
				errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',1.5, 'Color','k');
				line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 1.5);
				plot(data_MTF.fms, data_MTF.rate_sm, 'Color', 'b', 'LineWidth',2);
				hold off
				set(gca, 'XScale', 'log');
				xlim([data_MTF.fms(1) data_MTF.fms(end)])
				xticks([1 2 5 10 20 50 100 200 500])
				xlabel('Modulation Frequency')
				ylimit = ylim;
				ylim([0 ylimit(2)])
				grid on
				box on
				title(['MTF: ' data_MTF.MTF_shape ', 200Hz=' data_MTF.at_200])
				ylabel('Avg. Rate (sp/s)')
			end

			% Plot STRF
			nexttile
			STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
			imagesc(data_STRF.t, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
			set(gca,'Ydir','normal','XLim',data_STRF.tlims,'YLim',[param_ST.fpeaks(2) param_ST.fpeaks(end)]./1000)
			hold on
			yline(CF/1000, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			colormap(redblue)
			grid on
			xlabel('Time (s)');
			ylabel('Frequency (kHz)')

			% Plot real response
			nexttile
			hold on
			errorbar(data_ST.fpeaks,data_ST.rate,data_ST.rate_std/(sqrt(param_ST.nrep)), 'LineWidth',1.5);
			plot(data_ST.fpeaks,(STRFmodel.avModel.*STRFmodel.ratio), 'LineWidth',1.5);
			plot(data_ST.fpeaks, RM_R2.rate_RM, 'LineWidth',1.5);
			xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
			xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
			grid on
			xlabel('Tone Frequency (Hz)')
			ylabel('Avg Rate (sp/s)')
			ylim([0 STRFmodel.max_all+7])
			xticklabels(xticks/1000)
			title('STRF and RM Predictions')
			hLeg = legend({'Data', sprintf('STRF, R^2=%0.2f', STRFmodel.R2),...
				sprintf('RM, R^2=%0.2f',RM_R2.R2)}, 'Location','best', 'FontSize',5);
			hLeg.ItemTokenSize = [6,6];


			%% Add plots to PDF

			label = sprintf("%s, CF = %0.0fHz\n", putative, CF);
			p = Paragraph(label);
			p.FontSize = "14pt";
			p.WhiteSpace = "preserve";
			append(rpt,p);

			[plt1, images] = addtoSTPDF(images, fig, [putative '_' SPLs{ispl}]);
			append(rpt, plt1);

			label = '';
			p = Paragraph(label);
			p.FontSize = "14pt";
			p.WhiteSpace = "preserve";
			append(rpt,p);

		end
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
values = [8, 2.5];
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
