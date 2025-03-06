%% PDF_dog_v_gaussian_fits
clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load and initialize

% Load in spreadsheet
[base, datapath, savepath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "dog_analysis")


% Initialize report
filename = 'DoG_Test';
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


%%

% Sort by CF
num_sessions = length(dog_analysis);
linewidth = 0.75;
fontsize = 10;

% Plot each neuron
for isesh = 1:num_sessions

	% Load in data
	putative = dog_analysis(isesh).putative;
	CF = dog_analysis(isesh).CF;
	fpeaks = dog_analysis(isesh).fpeaks;

	% Paragraph intro
	label = sprintf("%s, CF = %0.0fHz\n", putative, CF);
	p = Paragraph(label);
	p.FontSize = "14pt";
	p.WhiteSpace = "preserve";
	append(rpt,p);

	% Output
	fprintf('Creating plots... %s, CF = %0.0fHz\n', putative, CF);

	% Plot data
	fig = figure;
	tiledlayout(1, 2)

	nexttile
	hold on
	errorbar(fpeaks./1000, dog_analysis(isesh).rate, dog_analysis(isesh).rate_std/sqrt(30), ...
		'linewidth', linewidth, 'color', 'b')
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	yline(dog_analysis(isesh).spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
	ylabel('Avg. Rate (sp/s)')
	xlabel('Spectral Peak Freq. (kHz)')

	% Plot gaussian
	plot(fpeaks./1000, dog_analysis(isesh).gaus_predicted, 'r', 'linewidth', linewidth)
	gaussian_adj_r_squared = dog_analysis(isesh).R2_gauss;

	% Plot DoG
	plot(fpeaks./1000, dog_analysis(isesh).dog_predicted, 'g', 'linewidth', linewidth)
	dog_adj_r_squared = dog_analysis(isesh).R2_dog;
	legend('Data', 'CF', 'Spont', 'Gaussian', 'DoG', 'location', 'westoutside')
	set(gca, 'FontSize',fontsize)

	% Comparing based on how close the curves are to data
	title(sprintf('Gaussian vs DoG Fits, p=%0.4f', dog_analysis(isesh).p_value))

	% Annotations
	gaus_msg = sprintf('Gaussian R^{2}=%0.02f', gaussian_adj_r_squared);
	text(0.05, 0.97, gaus_msg, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)
	dog_msg = sprintf('DoG R^{2}=%0.02f', dog_adj_r_squared);
	text(0.05, 0.89, dog_msg, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)

	% Plot DoG Parameters
	nexttile
	Fs = 100000;
	f = linspace(0, Fs/2, 100000);
	s_exc = dog_analysis(isesh).dog_params(1);
	s_inh = dog_analysis(isesh).dog_params(2);
	sigma_exc = 10^dog_analysis(isesh).dog_params(3);
	sigma_inh = 10^dog_analysis(isesh).dog_params(4);
	CF_exc = 10^dog_analysis(isesh).dog_params(5);
	CF_inh = 10^dog_analysis(isesh).dog_params(6);
	gauss_exc = normpdf(f, CF_exc, sigma_exc);
	gauss_inh = normpdf(f, CF_inh, sigma_inh);
	gauss_exc = s_exc*(gauss_exc./max(gauss_exc));
	gauss_inh = s_inh*(gauss_inh./max(gauss_inh));
	W = gauss_exc - gauss_inh;
	plot(gauss_exc)
	hold on
	plot(-1*gauss_inh)
	plot(W)
	title('DoG Filter')
	set(gca, 'fontsize', fontsize)
	ylabel('Amplitude')
	xlabel('Frequency (Hz)')
	xline(CF, '--', 'linewidth', 2)
	set(gca, 'xscale', 'log')
	xlim([300 10000])
	xticks([200 500 1000 2000 5000 10000])
	xticklabels([200 500 1000 2000 5000 10000]./1000)

	% Add to PDF
	[plt1, images] = addtoSTPDF(images, fig, putative);
	append(rpt, plt1);
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
values = [8, 3];
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
