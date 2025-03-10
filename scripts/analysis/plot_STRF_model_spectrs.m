%% strf_model_predictions
%
%
% J. Fritzinger, created 10/6/24

%% Load in spreadsheet

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Run model for all units

binmodes = {'Contra', 'Binaural'};
SPLs = {'63', '73'};
ispl = 1;
ibin = 2;

% Find sessions of interest
has_data = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
index = find(has_data);
isesh = 10;

for itype = 1:3

	% Load in data
	s_ind = index(isesh);
	putative_neuron = sessions.Putative_Units{s_ind};
	CF = sessions.CF(s_ind);
	load(fullfile(datapath, 'neural_data', [putative_neuron '.mat']), 'data');
	params_STRF = data{4,ibin};
	param_ST = data(6+ispl,ibin);

	% General analysis
	data_STRF = analyzeSTRF(params_STRF);
	data_ST = analyzeST(param_ST, CF);
	param_ST = param_ST{1};
	data_ST = data_ST{1};

	% Calculate model response
	param_ST.Fs = 48000;
	param_ST.mnrep = param_ST.nrep;
	param_ST.dur = 0.3;
	[R2, avModel, stdModel, ratio, max_all] = modelTimbreSTRF(param_ST, data_STRF, data_ST, itype);

	% Display progress
	fprintf('%s Done, %.2f through %s \n', putative_neuron, isesh/length(index), binmodes{ibin})

	%% Plot
	figure('Position',[3,632,623,273])
	tiledlayout(1, 2)

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
	plot(data_ST.fpeaks,(avModel.*ratio), 'LineWidth',1.5);
	xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
	xlim([param_ST.fpeaks(1) param_ST.fpeaks(end)]);
	grid on
	xlabel('Tone Frequency (Hz)')
	ylabel('Avg Rate (sp/s)')
	ylim([0 max_all+7])
	xticklabels(xticks/1000)
	title(sprintf('R^2 = %0.2f\n', R2))

end
