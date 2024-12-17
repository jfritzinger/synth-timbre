%% New Spreadsheet with all Spectral Centroid data 
% J. Fritzinger, updated 1/11/24
%
clear

%% Load in stimuli_sessions 

% Reads in spreadsheet & load data
if ismac
	%filename = '/Volumes/CarneyLab/Rabbit_data/session_table.xlsx';
	filename = '/Volumes/Rabbit_data/session_table.xlsx';
else
	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
end
sessions = readtable(filename, 'PreserveVariableNames',true);

%% Find sessions of interest 

rabbit = [24, 25, 26, 27, 28, 29, 30];
rating = {'Good', 'Excellent'};
interest_rabbit = ismember(sessions.Rabbit,rabbit);
interest_rating = ismember(sessions.Rating,rating);

stimulus = 'SPEC_slide_Spectral_Centroid';
names = sessions.Properties.VariableNames;
ind = find(strcmp(names,stimulus));
stim_columns = table2array(sessions(:,ind));
stim_columns(stim_columns > 1) = 1; % ignores Paul's system of using 1 and 2
interest_stim = sum(stim_columns,2);

interest = interest_rabbit & interest_rating & interest_stim;
index = find(interest);
num_interest = length(index);


%% Create empty table 

var_names = ["Session","TT","N", "CF", "MTF", "BMF", "WMF", "STRF", "Version", ...
	"43dB", "63dB", "73dB", "83dB", "43dB_100", "63dB_100", "73dB_100", "83dB_100",...
	"43dB_con", "63dB_con", "73dB_con", "83dB_con","Error"];
var_types = ["string","double","double", "double", "string", "double", "double", "double", "double", ...
	"double", "double", "double", "double", "double", "double", "double", "double", ...
	"double", "double", "double", "double", "string"];
num_vars = length(var_names);
listData = table('Size',[num_interest, num_vars], 'VariableTypes', var_types,...
    'VariableNames',var_names);

%% Fill out table (from post_process)

[userid, base_dir, ~, report_path, data_path] = findPaths();
for iclus = 1:num_interest

	% Fill out table from sessions 
	session = sprintf('R%03dS%03d', sessions.Rabbit(index(iclus)), sessions.Session(index(iclus)));
	listData.Session(iclus) = session;
	listData.TT(iclus) = sessions.Tetrode(index(iclus));
	listData.N(iclus) = sessions.Neuron(index(iclus));
	listData.CF(iclus) = sessions.CF(index(iclus));
	listData.MTF(iclus) = sessions.MTF(index(iclus));
	listData.BMF(iclus) = sessions.BMF(index(iclus));
	listData.WMF(iclus) = sessions.WMF(index(iclus));
	listData.STRF(iclus) = sessions.STRF(index(iclus));
	

	% Get the path to each session
	rab_num = num2str(sessions.Rabbit(index(iclus)));
	CF = sessions.CF(index(iclus));
	if ismac
		session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
	else
		session_dir_name = base_dir{contains(base_dir, rab_num)};
	end
	session_dir = fullfile(session_dir_name, session);

	% Load in session, or add message saying it could not be loaded
	try
		[clusters, params, ~] = loadPhysiologySession(session_dir, session, userid);
	catch
		msg = sprintf('Would not load in post_process');
		listData.Error(iclus) = msg;
		continue
	end
	cluster = clusters([clusters.tetrode] == sessions.Tetrode(index(iclus))); % Select tetrode
	cluster = cluster([cluster.neuron] == sessions.Neuron(index(iclus))); % Select neuron

	% Check that the version is 4+
	has_CF = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
		strcmp(p.SPEC_slide_type,'Spectral_Centroid'), params);
	params_st = params(has_CF);
	version = cellfun(@(p) p.version, params_st);
	version = unique(version);
	listData.Version(iclus) = version;
	if version <= 3
		msg = sprintf('Version %d', version);
		listData.Error(iclus) = msg;
		continue
	end

	% Match to CF
	if any(CF) 
		has_CF = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
			strcmp(p.SPEC_slide_type,'Spectral_Centroid'), params);
		if any(has_CF)
			fpeak_min = cellfun(@(p)p.fpeaks(1), params(has_CF));
			fpeak_max = cellfun(@(p)p.fpeaks(end), params(has_CF));
			range = [min(fpeak_min) max(fpeak_max)];
			if CF <= min(range,[],'all') || CF >= max(range,[],'all')
				msg = sprintf('Does not match sliding CF');
				listData.Error(iclus) = msg;
				continue
			end
		end
	end

	% Get data for rest of table 
	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, near_CF]...
		= finddata(params, CF);

	if ~any(F0_100 | F0_200)
		msg = 'Wrong F0';
		listData.Error(iclus) = msg;
		continue
	end

	bin_200 = [any(bin & F0_200 & level_43 & even) any(bin & F0_200 & level_63 & even)...
		any(bin & F0_200 & level_73 & even) any(bin & F0_200 & level_83 & even)];
	bin_100 = [any(bin & F0_100 & level_43 & even) any(bin & F0_100 & level_63 & even)...
		any(bin & F0_100 & level_73 & even) any(bin & F0_100 & level_83 & even)];
	contra_200 = [any(contra & F0_200 & level_43 & even) any(contra & F0_200 & level_63 & even)...
		any(contra & F0_200 & level_73 & even) any(contra & F0_200 & level_83 & even)];

	% Fill out table 
	for i = 1:4
		if i == 1 && bin_200(i) == 1
			listData.("43dB")(iclus) = any(bin & F0_200 & level_43 & even);
		elseif i == 2 && bin_200(i) == 1
			listData.("63dB")(iclus) = any(bin & F0_200 & level_63 & even);
		elseif i == 3 && bin_200(i) == 1
			listData.("73dB")(iclus) = any(bin & F0_200 & level_73 & even);
		elseif i == 4 && bin_200(i) == 1
			listData.("83dB")(iclus) = any(bin & F0_200 & level_83 & even);
		end
	end

	for i = 1:4
		if i == 1 && bin_100(i) == 1
			listData.("43dB_100")(iclus) = any(bin & F0_100 & level_43 & even);
		elseif i == 2 && bin_100(i) == 1
			listData.("63dB_100")(iclus) = any(bin & F0_100 & level_63 & even);
		elseif i == 3 && bin_100(i) == 1
			listData.("73dB_100")(iclus) = any(bin & F0_100 & level_73 & even);
		elseif i == 4 && bin_100(i) == 1
			listData.("83dB_100")(iclus) = any(bin & F0_100 & level_83 & even);
		end
	end

	for i = 1:4
		if i == 1 && contra_200(i) == 1
			listData.("43dB_con")(iclus) = any(contra & F0_200 & level_43 & even);
		elseif i == 2 && contra_200(i) == 1
			listData.("63dB_con")(iclus) = any(contra & F0_200 & level_63 & even);
		elseif i == 3 && contra_200(i) == 1
			listData.("73dB_con")(iclus) = any(contra & F0_200 & level_73 & even);
		elseif i == 4 && contra_200(i) == 1
			listData.("83dB_con")(iclus) = any(contra & F0_200 & level_83 & even);
		end
	end

end

%% Save spreadsheet

writetable(listData,'SynthSessions.xlsx')


%% Functions 

function [bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, near_CF] = finddata(population, CF)

	% Find binaural
	bin = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==2, population);

	% Find contra
	contra = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==1, population);

	% Find F0 = 100Hz
	F0_100 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==100, population);

	% Find F0 = 200Hz
	F0_200 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==200, population);

	% Find 43 db SPL
	level_43 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==43 || p.spl==40), population);

	% Find 63 dB SPL
	level_63 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==63 || p.spl==60), population);

	% Find 73 dB SPL
	level_73 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==73 || p.spl==70), population);

	% Find 83 dB SPL
	level_83 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==83 || p.spl==80), population);

	% Find even 200Hz fpeak_mid
	even = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);

	% Find nearest to CF 
	datasets = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);
	near_CF = cellfun(@(p) min(abs(p.fpeak_mid-CF)), population(datasets));
end