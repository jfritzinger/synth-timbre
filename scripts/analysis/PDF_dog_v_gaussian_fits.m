%% PDF_dog_v_gaussian_fits
%clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load and initialize
if ismac
	addpath('/Users/jfritzinger/Projects/shared-models/DoG-model', '-begin')
else
	addpath('C:\Projects_JBF\shared-models\DoG-model', '-begin')
end

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
filename = 'DoG_Gauss_Compareadd ';
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

has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);
linewidth = 0.75;
fontsize = 5;

% Plot each neuron
R2_dog_all = NaN(1, num_sessions);
R2_gauss_all = NaN(1, num_sessions);
R2_dog_all2 = NaN(1, num_sessions);
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

		% Calculate fits and plot

		% RM to get spont
		params_RM = data{2, 2};
		data_RM = analyzeRM(params_RM);
		spont = data_RM.spont;

		% Synthetic timbre analysis
		params = data(7, 2);
		params = params(~cellfun(@isempty, params));
		data_ST  = analyzeST(params, CF);
		data_ST = data_ST{1};
		rate = data_ST.rate;
		rate_std = data_ST.rate_std;
		rlb = data_ST.rlb;
		rub = data_ST.rub;
		fpeaks = data_ST.fpeaks;
		spl = data_ST.spl;
		rate_sm = data_ST.rates_sm;
		max_rate = max(rate);

		% Generate stimulus
		params{1}.Fs = 100000;
		params{1}.physio = 1;
		params{1}.mnrep = 1;
		params{1}.dur = 0.3;
		params{1} = generate_ST(params{1});
		params{1}.num_stim = size(params{1}.stim, 1);

		% Run models 
		Fs = 100000;
		observed_rate = rate;
		r0 = spont;
		stim = params{1}.stim;
		%[gaussian_params, dog_params, dog_params2] = fitGaussAndDoG(params, CF, Fs, observed_rate, r0);
		nrep = 15;
		num_params = 5;
		dog_params = fit_dog_model(nrep, CF, Fs, stim, observed_rate, r0, num_params);
		nrep = 5;
		num_params = 6;
		dog_params6 = fit_dog_model_6param(nrep, CF, Fs, stim, observed_rate, r0);
		nrep = 10;
		gaussian_params = fit_gaussian_model(nrep, CF, Fs, stim, observed_rate, r0);

		% Plot data
		fig = figure;
		tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact')
		nexttile
		hold on
		errorbar(fpeaks./1000, rate, rate_std/sqrt(params{1}.nrep), ...
			'linewidth', linewidth, 'color', 'k')
		xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
		yline(spont, 'color', [0.5 0.5 0.5], LineWidth=linewidth)
		ylabel('Avg. Rate (sp/s)')
		xlabel('Spectral Peak Freq. (kHz)')

		% Plot gaussian
		stim = params{1}.stim;
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		gaus_predicted = zeros(nstim, 1);
		for i = 1:nstim
			W = gaussian_model(f, gaussian_params);
			gaus_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end
		plot(fpeaks./1000, gaus_predicted, 'r', 'linewidth', linewidth)
		gaussian_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			gaus_predicted, 3);
		gaus_msg = sprintf('Gaussian,\n R^2=%0.02f', gaussian_adj_r_squared);

		% Plot DoG
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		dog_predicted = zeros(nstim, 1);
		for i = 1:nstim
			W = dog_model(f, dog_params);
			dog_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end
		plot(fpeaks./1000, dog_predicted, 'g', 'linewidth', linewidth)
		dog_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			dog_predicted, 5);
		dog_msg = sprintf('DoG,\n R^2=%0.02f', dog_adj_r_squared);

		% Plot DoG2
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		dog_predicted2 = zeros(nstim, 1);
		for i = 1:nstim
			W = dog_model(f, dog_params6);
			dog_predicted2(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end
		plot(fpeaks./1000, dog_predicted2, 'b', 'linewidth', linewidth)
		dog_adj_r_squared2 = calculate_adj_r_squared(observed_rate,...
			dog_predicted2, 6);
		set(gca, 'FontSize',fontsize)
		dog_msg2 = sprintf('DoG6,\n R^2=%0.02f', dog_adj_r_squared2);

		hleg = legend('', '', '', gaus_msg, dog_msg, dog_msg2, 'location', 'westoutside');
		hleg.ItemTokenSize = [4, 4];

		% Comparing based on how close the curves are to data
		p_value = ftest(rate, gaus_predicted, dog_predicted);
		title(sprintf('Gaussian vs DoG Fits, p=%0.4f', p_value))
		set(gca, 'FontSize',fontsize)

		% Plot DoG Parameters
		nexttile
		hold on

		Fs = 100000;
		f = linspace(0, Fs/2, 100000);
		W = gaussian_model(f, gaussian_params);
		plot(f/1000,W, 'color', 'r')

		W2 = dog_model(f, dog_params);
		plot(f/1000,W2, 'color', 'g')

		W3 = dog_model(f, dog_params6);
		plot(f/1000,W3, 'color', 'b')

		% Plot labels
		xline(CF/1000, '--', 'linewidth', 1.5)
		title('DoG and Gaussian Kernels')
		set(gca, 'fontsize', fontsize)
		ylabel('Amplitude')
		xlabel('Frequency (kHz)')
		xlim([200 10000]/1000)
		set(gca, 'xscale', 'log')
		grid on
		xticks([0.1 0.2 0.5 1 2 5 10])

		% Get R^2 for all
		R2_dog_all(isesh) = dog_adj_r_squared;
		R2_gauss_all(isesh) = gaussian_adj_r_squared;
		R2_dog_all2(isesh) = dog_adj_r_squared2;

		% Get f-test for all
		p_value = ftest(rate, gaus_predicted, dog_predicted);
		p_value2 = ftest(rate, gaus_predicted, dog_predicted2);

		% Struct to save out all data and fits 
		dog_gauss_analysis.putative = putative;
		dog_gauss_analysis.dog_predicted = dog_predicted;
		dog_gauss_analysis.dog_predicted2 = dog_predicted2;
		dog_gauss_analysis.gaus_predicted = gaus_predicted;
		dog_gauss_analysis.CF = CF;
		dog_gauss_analysis.rate = observed_rate;
		dog_gauss_analysis.R2_dog = dog_adj_r_squared;
		dog_gauss_analysis.R2_dog2 = dog_adj_r_squared2;
		dog_gauss_analysis.R2_gauss = gaussian_adj_r_squared;
		dog_gauss_analysis.fpeaks = data_ST.fpeaks;
		dog_gauss_analysis.spont = spont;
		dog_gauss_analysis.rate_std = data_ST.rate_std;
		dog_gauss_analysis.p_value_gd = p_value;
		dog_gauss_analysis.p_value_gd2 = p_value;
		dog_gauss_analysis.dog_params = dog_params;
		dog_gauss_analysis.dog_params2 = dog_params6;
		dog_gauss_analysis.gauss_params = gaussian_params;

		filename = [putative '.mat'];
		%savepath = '/Volumes/Synth-Timbre/data/manuscript/';
		savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript';
        %savepath = '\\NSC-LCARNEY-H2\DataFiles_JBF\Synth-Timbre\data\manuscript';
		save(fullfile(savepath, 'dog_model', filename), 'dog_gauss_analysis')

		% Struct to save out all data and fits 
		dog_analysis(isesh).putative = putative;
		dog_analysis(isesh).dog_predicted = dog_predicted;
		dog_analysis(isesh).dog_predicted2 = dog_predicted2;
		dog_analysis(isesh).gaus_predicted = gaus_predicted;
		dog_analysis(isesh).CF = CF;
		dog_analysis(isesh).rate = observed_rate;
		dog_analysis(isesh).R2_dog = dog_adj_r_squared;
		dog_analysis(isesh).R2_dog2 = dog_adj_r_squared2;
		dog_analysis(isesh).R2_gauss = gaussian_adj_r_squared;
		dog_analysis(isesh).fpeaks = data_ST.fpeaks;
		dog_analysis(isesh).spont = spont;
		dog_analysis(isesh).rate_std = data_ST.rate_std;
		dog_analysis(isesh).p_value = p_value;
		dog_analysis(isesh).dog_params = dog_params;
		dog_analysis(isesh).dog_params2 = dog_params6;
		dog_analysis(isesh).gauss_params = gaussian_params;

		fprintf('%s done, %d percent done\n', putative, round(isesh/num_sessions*100))

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

%% 

save(fullfile(datapath, 'dog_analysis3.mat'), "dog_analysis", ...
	"R2_gauss_all", "R2_dog_all", "R2_dog_all2")

%% FUNCTIONS

function [img, images] = addtoSTPDF(images, fig, title)
import mlreportgen.dom.*

% Set figure size, recommended
values = [4.2, 1.5];
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
