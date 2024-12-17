function pdf_grid_search(putative_timbre, CF)

% Set up save paths
[~, computer] = system('hostname');
if ismac
	savepath = '/Volumes/Synth-Timbre/data/manuscript/model-fits';
	addpath('/Users/jfritzinger/Projects/shared-models/efferent-model/')
elseif contains(computer, 'I1') % I1
	savepath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
else
	savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript\model-fits';
end

% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
report_name = fullfile(savepath, putative_timbre, ...
	sprintf('%s_Model.pdf', putative_timbre));
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.05in';
pm.PageMargins.Header = '0.05in';
pm.PageMargins.Bottom = '0.05in';
pm.PageMargins.Footer = '0.05in';
pm.PageMargins.Left = '0.05in';
pm.PageMargins.Right = '0.05in';

num_paramCF = 5;
for iparamCF = 1:num_paramCF

	% Load in IC model
	filename = sprintf('%s_IC_%d.mat', putative_timbre, iparamCF);
	load(fullfile(savepath, putative_timbre, filename), ...
		'params', 'model_params', 'lateral_model', 'paramS_range',  ...
		'paramCF_range')
	num_paramS = length(paramS_range);

	% Generate page heading
	title_pdf = sprintf('%s, CF range: %0.2g', ...
		putative_timbre, paramCF_range(iparamCF));
	h = Heading(2, title_pdf);
	b = Border();
	b.BottomStyle = 'single';
	b.BottomColor = 'LightGray';
	b.BottomWidth = '1pt';
	h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	append(rpt,h);

	% Plot data
	fig = figure();
	tiledlayout(1, 4, 'TileSpacing','compact', 'Padding','tight')
	linewidth = 1.5;

	% Plot MTF
	data_MTF = analyzeMTF(params{1});
	nexttile()
	hold on
	line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',...
		[0.4 0.4 0.4], 'LineWidth',linewidth);
	errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params{1}.nrep),...
		'.k', 'LineWidth',linewidth, 'CapSize',4);
	line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', ...
		'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
	hold off
	set(gca, 'XScale', 'log');
	xlim([data_MTF.fms(1) data_MTF.fms(end)])
	xticks([1 2 5 10 20 50 100 200 500])
	xlabel('Modulation Frequency (Hz)')
	grid on
	box on
	ylabel('Avg. Rate (sp/s)')
	title('MTF')
	set(gca, 'fontsize', 6)

	% Plot WB-TIN
	for iWB = 1:3
		if ~isempty(params{1+iWB})
			% Analysis
			data_WB = analyzeWBTIN(params(iWB+1), []);
			data_WB = data_WB{1};

			% Plot
			nexttile
			spont_color = [0.4 0.4 0.4];
			colors = {'#82BB95', '#3F985C', '#034E1C'};
			hold on
			yline(mean(data_WB.rate(1))/1000,'Color',[0.4 0.4 0.4], 'LineWidth', 1.5);
			xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.5); % CF line
			errorbar(data_WB.fpeaks(:,1)/1000,data_WB.rate(:,1),...
				data_WB.rate_std(:,1)/(sqrt(params{iWB+1}.nrep)), ...
				'LineWidth', 2,'Color', colors{2});
			errorbar(data_WB.fpeaks(:,2)/1000,data_WB.rate(:,2),...
				data_WB.rate_std(:,2)/(sqrt(params{iWB+1}.nrep)), ...
				'LineWidth', 2,'Color', colors{3});
			plot_range = [params{iWB+1}.fpeaks(2) params{iWB+1}.fpeaks(end)];
			set(gca, 'XScale', 'log');
			box on
			grid on
			xlabel('Tone Frequency (Hz)')
			ylabel('Avg Rate (sp/s)')
			title(['No = ' num2str(params{iWB+1}.No)])
			set(gca, 'fontsize', 6)
		else
			nexttile
			title('No data')
			set(gca, 'fontsize', 6)
		end
	end

	% Append
	tempname = sprintf('%s_CF%d_data', putative_timbre, iparamCF);
	[plt, images] = addToPDF(images, fig, tempname, [8, 1.25]);
	append(rpt, plt);

	for iparamS = 1:num_paramS

		% Add paragraph of stimulus parameters and outputs
		param_out = sprintf('CF range = %0.3g, Strength = %0.3g, Delay = 0',...
			paramCF_range(iparamCF), paramS_range(iparamS));
		p = Paragraph(param_out);
		append(rpt,p);

		fig = figure();
		tiledlayout(1, 4, 'TileSpacing','compact', 'Padding','tight')

		% Plot MTF
		nexttile
		spont_color = [0.4 0.4 0.4];
		[~, avIC, stdIC] = plotMTF(params{1}, lateral_model{iparamS, 1}.avIC, 0);
		yline(avIC(1), 'Color',spont_color, 'LineWidth',2)
		hold on
		errorbar(params{1}.all_fms, avIC, stdIC./sqrt(params{1}.mnrep), 'LineWidth',2,'Color','k')
		y = ylim;
		ylim([0,y(2)+5])
		xlim([params{1}.all_fms(2) params{1}.all_fms(end)])
		set(gca,'xtick',[1.2,2, 5,  20, 50, 200, 500],'xticklabel',...
			{'Unmod','2','5','20', '50','200','500'},'xscale','log')
		set(gca, 'XScale', 'log');
		ylabel('Rate (sp/s)','fontsize',10)
		xlabel('Modulation Frequency (Hz)')
		box on
		grid on
		set(gca, 'fontsize', 6)

		for iWB = 1:3 % Plot model
			if ~isempty(params{1+iWB})
				[SNRs,~,si] = unique([params{1+iWB}.mlist.SNR].');
				num_SNRs = length(SNRs);
				[fpeaks,~,fi] = unique([params{1+iWB}.mlist.fpeak].');
				num_fpeaks = length(fpeaks);
				rate_size = [num_fpeaks,num_SNRs];
				[avIC,stdIC,~,~] = accumstats({fi,si},lateral_model{iparamS, 1+iWB}.avIC, rate_size);

				nexttile
				yline(avIC(1,1), 'k', 'LineWidth',1.5)
				hold on
				errorbar(params{1+iWB}.fpeaks/1000, avIC, stdIC./sqrt(params{1+iWB}.mnrep),...
					'LineWidth',1.5)
				xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.5); % CF line
				xlim([params{1+iWB}.freq_lo params{1+iWB}.freq_hi]./1000)
				set(gca, 'XScale', 'log');
				xlabel('Tone Frequency (kHz)')
				box on
				grid on
				y = ylim;
				ylim([0,y(2)])
				set(gca, 'fontsize', 6)
			else
				nexttile
				set(gca, 'fontsize', 6)
			end
		end

		% Append
		tempname = sprintf('%s_CF%d_param%d', putative_timbre, iparamCF, iparamS);
		[plt, images] = addToPDF(images, fig, tempname, [8, 1.25]);
		append(rpt, plt);

		% Show progress 
		fprintf('Model CF range %d, %d/%d done \n', iparamCF, iparamS, num_paramS)

	end
end

close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)

end

%% 

function [img, images] = addToPDF(images, fig, title, size)
import mlreportgen.dom.*

% Set figure size, recommended
fig.PaperSize = size;
fig.PaperPosition = [0 0 size];
fig.Units = 'inches';
fig.Position(3:4) = size;

% Add the plot to the document
if ismac
	temppath = '/Users/jfritzinger/Documents/MATLAB/Temp_PDF';
else
	temppath = 'C:\Users\jfritzinger\Documents\MATLAB\Temp_PDF';
end
name = sprintf('%s.svg', title);
tempfilepath = fullfile(temppath, name);
print(fig, tempfilepath, '-dsvg');
img = Image(tempfilepath);
delete(fig) %delete plot figure window
images = [images {img}];

end