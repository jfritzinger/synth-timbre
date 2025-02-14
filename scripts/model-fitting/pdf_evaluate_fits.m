function pdf_evaluate_fits(putative, CF, putative_timbre)

%% evaluate_fits

[~, computer] = system('hostname');
if ismac
	modelpath = '/Volumes/WB-TIN/data/model-fits';
	addpath('/Users/jfritzinger/Projects/WB-TIN/scripts/helper-functions',...
		'-end')
elseif contains(computer, 'I1') % I1
    modelpath = '\\NSC-LCARNEY-H2\Synth-Timbre\data\manuscript\model-fits';
else
	modelpath = 'C:\DataFiles_JBF\WB-TIN\data\model-fits';
	addpath('C:\Projects_JBF\WB-TIN\scripts\helper-functions\', '-end')
end
[~, datapath, ~, ~] = getPathsWBTIN();

%% Load in data 

load(fullfile(datapath, [putative '.mat']), 'data');
data_rates = analyze_data(data, CF); % Analyze data and put in correct form 

%% Load in AN response

% Load in AN
filename = sprintf('%s_AN.mat', putative_timbre);
load(fullfile(modelpath, putative_timbre, filename), 'params', 'AN', 'model_params')

% Load in IC parameter values 
filename = sprintf('%s_IC.mat', putative_timbre);
load(fullfile(modelpath, putative_timbre, filename), 'fit_params_all')

%% Set up PDF

% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
if ismac
	report_path = '/Volumes/WBTIN_ModelFits';
else
	report_path = 'C:\DataFiles_JBF\WBTIN_ModelFits';
end

filename = 'Pearson';
report_name = fullfile(modelpath, putative_timbre, ...
	sprintf('%s_%s.pdf', putative_timbre, filename));
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

% Evaluate model fits  
data_MTF = data_rates(1:26);
data_WB = data_rates(27:end);

num_paramCF = size(AN, 1);
for iparamCF = 1:num_paramCF

	% Run model with fit parameters
	AN_sub = AN(iparamCF,:);
	fit_params = fit_params_all(iparamCF,:);
	CS_params = [fit_params(1:2) 0.001];
	BMFs = fit_params(3:5);
	nstim = size(params, 2);
	model_outputs = cell(nstim, 1);
	for istim = 1:nstim
		param = params{istim};
		an_sout = squeeze(AN_sub{istim}.an_sout);
		an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
		an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
		model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
			an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
			'BMFs', BMFs);
	end

	% Analyze model output
	% iWB = 1;
	for istim = 1:nstim
		param = params{istim};
		model_output = model_outputs{istim};
		switch param.type
			case 'typMTFN'
				[~, model_MTF, ~, ~] = plotMTF(param, model_output.avIC, 0);
			case 'TIN'
				[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
			case 'SPEC_slide' % WB-TIN only for now
				[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
				%model_WB = model_WB(:,2);
		end
	end

	% Calculate MSE 
	mse = minimize_IC_model_fit(data_rates, AN_sub, params, model_params,...
		fit_params);
	mse_all(iparamCF) = mse;

	% Paragraph 
	msg = sprintf(['CF range: %0.2f, [%0.02f %0.02f 0.001], ...' ...
		'BMF=[%0.0f %0.0f %0.0f], R=%0.2f'], ...
		AN_sub{1,1}.CF_span, fit_params(1), fit_params(2), fit_params(3),...
		fit_params(4), fit_params(5), mse);
	p = Paragraph(msg);
	p.FontSize = "14pt";
	p.WhiteSpace = "preserve";
	append(rpt,p);

	fig = figure();
	tiledlayout(1, 3, 'Padding','compact')

	% Plot MTF
	nexttile
	param = params{1};
	hold on
	yline(data_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,data_MTF);
	yline(model_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,model_MTF);
	%plot(data.fms,rate_sm,'-b', 'LineWidth', 1)
	hold off
	xtick = [1 2 5 20 50 100 200 500];
	xlim(xtick([1 end]))
	xlabel('Modulation Freq (Hz)')
	ylabel('Avg. Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	legend('Unmodulated', 'Location','best')
	grid on
	axis_set = axis;
	axis_set(3) = 0;
	axis(axis_set);
	r = corrcoef(model_MTF, data_MTF);
	message = sprintf('R=%0.2f', r(1,2));
	text(0.05, 0.95, message, 'Units', 'normalized', ...
			'VerticalAlignment', 'top')

	% Plot NB-TIN
	nexttile
% 	hold on
% 	num_SNRs = 3;
% 	SNRs = [-inf params{2}.SNR];
% 	x = repmat(1:num_SNRs, 1, 1);
% 	plot(x, data_TIN, 'LineWidth',1.5);
% 	plot(x,model_TIN, 'LineWidth',1.5);
% 	xlabel('SNR');
% 	xticks(1:num_SNRs)
% 	xticklabels([-inf 30 40])
% 	xlim([0 num_SNRs+1])
% 	ylabel('Spike rate (sp/s)')
% 	y = ylim;
% 	ylim([0,y(2)])
% 	legend('Data', 'Model', 'Location','best')
% 	r = corrcoef(model_TIN, data_TIN);
% 	message = sprintf('R=%0.2f', r(1,2));
% 	text(0.05, 0.95, message, 'Units', 'normalized', ...
% 			'VerticalAlignment', 'top')

	nexttile
	param = params{2};
	freqs = param.fpeaks;
	hold on
	yline(mean(data_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,data_WB,  'LineWidth',2, 'Color',"#0072BD")
	yline(mean(model_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,model_WB,  'LineWidth',2, 'Color',"#D95319");
	xline(CF/1000, ':', 'Color', [0.5 0.5 0.5]);
	hold off
	xlabel('Tone Frequency (kHz)')
	ylabel('Avg Rate (sp/s)')
	set(gca, 'XScale', 'log');
	grid on
	box on
	xlim([param.fpeaks(2) param.fpeaks(end)]/1000);
	axis_set = axis;
	axis_set(3) = 0;
	axis(axis_set);
	r = corrcoef(model_WB, data_WB);
	message = sprintf('R=%0.2f', r(1,2));
	text(0.05, 0.95, message, 'Units', 'normalized', ...
		'VerticalAlignment', 'top')

	tempname = ['CF' num2str(iparamCF)];
	[plt, images] = addToPDF_model(images, fig, tempname, [8, 2.1]);
	append(rpt, plt);

	p = Paragraph('  ');
	p.FontSize = "14pt";
	p.WhiteSpace = "preserve";
	append(rpt,p);

end

% Close report
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)
end

