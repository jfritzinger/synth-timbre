%% save_peak_finding_results.m
%
% Script that...
%
%
% Author: J. Fritzinger
% Created: 2022-09-26; Last revision: 2024-09-26
%
% -------------------------------------------------------------------------
clear

%%

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Create table
varNames = ["Putative", "CF", "MTF","MTF_at200", "MTF_str", ...
	"SPL", "binmode", "F0", ...
	"Type", "Prom", "Width", "Lim", "Freq", "Q", "Q_log"];
varTypes = ["string", "double", "string", "string", "double", ...
	"double", "double", "double", ...
	"string", "double", "double", "double", "double", "double", "double"];
est_num_rows = 830; % set to number larger than
num_cols = length(varNames);
table_size = [est_num_rows num_cols];
tables = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

% Create has_data
data_ind(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
data_ind(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
data_ind(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
data_ind(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
data_ind(:,5) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_con);
data_ind(:,6) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_con);
data_ind(:,7) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB_con);
data_ind(:,8) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_con);
data_ind(:,9) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB_100);
data_ind(:,10) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB_100);
data_ind(:,11) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB_100);
has_data = any(data_ind,2);
data_index = find(has_data);
num_neurons = sum(has_data);

%% Plot each neuron
ii = 1;
for isesh = 1:num_neurons
	ineuron = data_index(isesh);

	% Load in session
	putative = sessions.Putative_Units{ineuron};
	CF = sessions.CF(ineuron);
	MTF_shape = sessions.MTF{ineuron};
	at200 = sessions.MTF_at200{ineuron};
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	if CF<2000
		CF_Group = 'Low';
	elseif CF>=2000 && CF<4000
		CF_Group = 'Med';
	else
		CF_Group = 'High';
	end

	% MTF analysis
	params_MTF = data{3, 2};
	if ~isempty(params_MTF)
		data_MTF = analyzeMTF(params_MTF);
	end

	for idata = 1:11
		if data_ind(ineuron,idata)==1

			% Load in proper dataset for each idata
			if ismember(idata, [1, 2, 3, 4])
				param_ST = data(5+idata, 2);
			elseif ismember(idata, [5, 6, 7, 8])
				param_ST = data(1+idata, 1);
			else
				param_ST = data(1+idata, 2);
			end

			% Analyze synthetic timbre
			spl = param_ST{1}.spl;
			data_ST = analyzeST(param_ST, CF);
			data_ST = data_ST{1};
			[~, peak_ind] = max(data_ST.rate);
			peak_f = data_ST.fpeaks(peak_ind);
			[~, peaksm_ind] = max(data_ST.rates_sm);
			peak_fsm = data_ST.fpeaks(peaksm_ind);

			% Find peaks & prominence values
			[peaks, dips, type, prom, width, lim, ~, ~, freq] = peakFinding(data_ST, CF);

			% Add data to table
			tables.Putative{ii} = sessions.Putative_Units{ineuron};
			tables.CF(ii) = CF;
			tables.CF_Group{ii} = CF_Group;
			tables.MTF{ii} = MTF_shape;
			tables.MTF_at200{ii} = at200;
			tables.MTF_str(ii) = data_MTF.perc_change;
			tables.SPL(ii) = spl;
			tables.Type{ii} = type;
			tables.binmode(ii) = param_ST{1}.binmode;
			tables.F0(ii) = param_ST{1}.Delta_F;
			tables.Width(ii) = width;
			tables.Lim(ii) = lim;
			tables.Prom(ii) = prom;
			tables.Freq(ii) = freq;
			tables.Q(ii) = freq/width;
			tables.Q_log(ii) = log10(freq/width);
			ii = ii+1;
		end
	end
	fprintf('%s done, %d percent done\n', putative, round(isesh/num_neurons*100))

end

%% Save table

% Save table
writetable(tables,'peak_picking2.xlsx')


