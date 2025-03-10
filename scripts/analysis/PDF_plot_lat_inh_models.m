%% plot_models_pdf.m
%
% Plots each synthetic timbre response from each neuron with the model
% responses for SFIE, energy, and population SFIE overlayed. 
%
%
% Author: J. Fritzinger
% Created: 2022-09-13; Last revision: 2024-09-13
%
% -------------------------------------------------------------------------
clear

import mlreportgen.dom.*
import mlreportgen.report.*

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Synth-Timbre/data/manuscript';
%modelpath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript';
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

% Initialize report
report_path = 'figures/pdfs/';
filename = 'Model_test';
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

%% Plot each dataset 

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
WB_TIN = cellfun(@(s)contains(s, 'R'), sessions.WBTIN_Units);
isMTF = strcmp(sessions.MTF, 'BE') | strcmp(sessions.MTF, 'BS');

bin200 = bin200 & WB_TIN & isMTF;
has_data = bin200(:,1) | bin200(:,2) | bin200(:,3) | bin200(:,4);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Plot each neuron
dips_p = NaN(num_sessions, 4);
peaks_p = NaN(num_sessions, 4);
for isesh = 1:40 %:num_sessions
	if isesh == 2 || isesh == 10 || isesh == 11 || isesh == 17
		continue
	end
	%ineuron = index(order(isesh));
	ineuron = index(isesh);
	if any(has_data(ineuron))

		% Load in data 
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))

		% Load in model data 
		load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
		lat_inh_orig = lat_inh;
		load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh2.mat']), 'lat_inh')

		% Paragraph intro
		label = '';
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		label = sprintf("%s, CF = %0.0fHz, %s\n", putative, CF, MTF_shape);
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		% Output 
		 fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);

		% Set up figure
		fig = figure('Position',[1,571,8*ppi,5*ppi]);
		tiledlayout(2, 4, 'TileSpacing','compact', 'Padding','tight')
		x_label = [1000 2000 3000 4000 6000 8000]./1000;
		fontsize = 7;

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
		xlabel('Modulation Frequency')
		set(gca,'fontsize',fontsize)
		ylimit = ylim;
		ylim([0 ylimit(2)])
		grid on
		box on
		title(['MTF: ' data_MTF.MTF_shape ', 200Hz=' data_MTF.at_200], 'FontSize', fontsize)
		ylabel('Avg. Rate (sp/s)')

		% Plot synthetic timbre (raw)
		spls = [43, 63, 73, 83];
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		for ispl = 1:4
			nexttile(4+ispl)
			if bin200(ineuron, ispl)==1
				param_ST = data(5+ispl, 2);
				data_ST = analyzeST(param_ST, CF);
				data_ST = data_ST{1};
				max_rate = max(data_ST.rate);
				spont = data_RM.spont;
				rate_subtracted = data_ST.rate - spont;

				% Z-score 
				rate = data_ST.rate; %zscore(data_ST.rate);
				rate_sm = data_ST.rates_sm; %zscore(data_ST.rates_sm);
				hold on
				plot(data_ST.fpeaks,rate, 'linewidth', 1.5, 'Color','k');
				yline(data_RM.spont, 'k')
				plot(data_ST.fpeaks, lat_inh_orig{ispl}.rate, 'LineWidth',1.5, 'Color','#CC79A7')
				plot(data_ST.fpeaks, lat_inh{ispl}.rate, 'LineWidth',1.5, 'Color','b')

				% Annotate lateral inhibition model R^2
				lefts = linspace(0.03, 0.8, 4);
				message = sprintf('R^2 lat inh = %.02f', lat_inh{ispl}.R2);
				annotation('textbox',[lefts(ispl) 0.05 0.2 0.0869],...
					'Color','k',...
					'String',message, ...
					'FontSize',8,...
					'EdgeColor','none');
				message = sprintf('R^2 lat inh orig = %.02f', lat_inh_orig{ispl}.R2);
				annotation('textbox',[lefts(ispl) 0.02 0.2 0.0869],...
					'Color','k',...
					'String',message, ...
					'FontSize',8,...
					'EdgeColor','none');

				% Calculate
				lat_inh_temp =  lat_inh{ispl}.rate;
				lat_inh_orig_temp = lat_inh_orig{ispl}.rate;

				% Normalize to 1 (will not use this when models are fit to rate correctly)
				lat_inh_temp = lat_inh_temp ./ max(lat_inh_temp);
				lat_inh_orig_temp = lat_inh_orig_temp ./ max(lat_inh_orig_temp);
				rate = data_ST.rate ./ max(data_ST.rate);
				p_l_l2 = ftest(rate, lat_inh_orig_temp, lat_inh_temp); % F-test LatInh/energy

				% Annotate 
				message = sprintf('p_l_l2=%.02f', p_l_l2);
				annotation('textbox',[lefts(ispl)+0.12 0.38 0.2 0.0869],...
					'Color','k',...
					'String',message, ...
					'FontSize',8,...
					'EdgeColor','none');
				plottitle = [num2str(spls(ispl)) ' dB SPL'];
				plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
				xlim(plot_range)
			else
				plottitle = [num2str(spls(ispl)) ' dB SPL'];
			end

			xline(CF, '--', 'Color',CF_color, 'linewidth', 1.5)
			xlabel('Frequency (Hz)')
			if ispl == 1
				ylabel('Z-score')
			end
			box on
			title(plottitle)
			grid on
			set(gca, 'fontsize', fontsize)
		end

		label = '';
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

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
values = [8, 4.65];
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
