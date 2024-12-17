%% New Spreadsheet with all Spectral Centroid data 
% J. Fritzinger, updated 1/11/24
clear

%% Load in stimuli_sessions 

% Reads in spreadsheet & load data
if ismac
	filename = '/Volumes/Rabbit_data/session_table.xlsx';
else
	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
end
sessions = readtable(filename, 'PreserveVariableNames',true);

%% Find sessions of interest 

rabbit = [27, 29];
interest_R24 = sessions.Rabbit==24 & sessions.Session>= 365;
interest_R25 = sessions.Rabbit==25 & sessions.Session>=356 & sessions.Session<631;
rating = {'Good', 'Excellent'};
interest_rabbit = ismember(sessions.Rabbit,rabbit);
interest_rating = ismember(sessions.Rating,rating);
driven = sessions.CF~=0;

stimulus = 'SPEC_slide_Spectral_Centroid';
names = sessions.Properties.VariableNames;
ind = find(strcmp(names,stimulus));
stim_columns = table2array(sessions(:,ind));
interest_stim = sum(stim_columns,2);

stimulus2 = 'Natural_Timbre';
ind = find(strcmp(names,stimulus2));
stim_columns2 = table2array(sessions(:,ind));
interest_stim2 = sum(stim_columns2,2);

interest = (interest_R24 | interest_R25 | interest_rabbit) & ...
	interest_rating & (interest_stim | interest_stim2) & driven;
index = find(interest);
num_interest = length(index);


%% Create empty table 

varNames = ["Rabbit", "Session","TT","N", "Putative_Units", ...
	"Include_SC", "Include_NT", "Pass", "Depth", "CF", "MTF", ...
	"BMF", "WMF", "MTF_con", "BMF_con", "WMF_con", "Error"];
varTypes = ["double", "double","double","double", "string", ...
	"string","string", "double", "double", ...
	"double", "string", "double", "double",...
	"string", "double", "double", "string"];
stimNames = ["char_spl", "char_ITD", "char_ILD", "type_RM", "type_RM_con", ...
	"typMTFN","typMTFN_con","STRF","STRF_con","SCHR" ...
	"43dB", "63dB", "73dB", "83dB", "43dB_100", "63dB_100", "73dB_100", "83dB_100",...
	"43dB_con", "63dB_con", "73dB_con", "83dB_con","Oboe", "Bassoon", "Other"];

stimTypes = repmat("double", 1, length(stimNames));
est_num_rows = 1000; % set to number larger than
num_cols = length([varNames stimNames]);
table_size = [est_num_rows num_cols];
listData = table('Size',table_size,'VariableTypes',[varTypes stimTypes],'VariableNames',[varNames stimNames]);


%% Fill out table (from post_process)

[userid, base_dir, ~, report_path, data_path] = findPaths();
for iclus = 1:num_interest

	% Fill out table from sessions 
	session = sprintf('R%03dS%03d', sessions.Rabbit(index(iclus)), sessions.Session(index(iclus)));
	
	listData.Rabbit(iclus) = sessions.Rabbit(index(iclus));
	listData.Session(iclus) = sessions.Session(index(iclus));
	listData.TT(iclus) = sessions.Tetrode(index(iclus));
	listData.N(iclus) = sessions.Neuron(index(iclus));
	listData.Pass(iclus) = sessions.Pass(index(iclus));
	listData.Depth(iclus) = sessions.Tetrode_Depth(index(iclus));

	listData.CF(iclus) = sessions.CF(index(iclus));
	listData.MTF(iclus) = sessions.MTF(index(iclus));
	listData.BMF(iclus) = sessions.BMF(index(iclus));
	listData.WMF(iclus) = sessions.WMF(index(iclus));
	listData.STRF(iclus) = sessions.STRF(index(iclus));
	listData.char_spl(iclus) = sessions.char_spl(index(iclus));
	listData.char_ITD(iclus) = sessions.char_ITD(index(iclus));
	listData.char_ILD(iclus) = sessions.char_ILD(index(iclus));
	listData.type_RM(iclus) = sessions.type_RM(index(iclus));
	listData.typMTFN(iclus) = sessions.typMTFN(index(iclus));

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
	cluster = clusters([clusters.tetrode] == sessions.Tetrode(index(iclus))); % Select tetrode
	cluster = cluster([cluster.neuron] == sessions.Neuron(index(iclus))); % Select neuron
	catch
		continue
	end

	 % Contra 
	 RM_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'type=RM') &&...
		 p.binmode==1, params);
	listData.type_RM_con(iclus) = any(RM_con);
	MTFN_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'typMTFN') &&...
		 p.binmode==1, params);
	listData.typMTFN_con(iclus) = any(MTFN_con);
	STRF_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'STRF') &&...
		 p.binmode==1, params);
	listData.STRF_con(iclus) = any(STRF_con);

	% Check that the version is 4+
	has_CF = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
		strcmp(p.SPEC_slide_type,'Spectral_Centroid'), params);	
	params_st = params(has_CF);
	if ~isempty(params_st)
		version = cellfun(@(p) p.version, params_st);
		version = unique(version);
		%listData.Version(iclus) = version;
		if version <= 3
			msg = sprintf('Version %d', version);
			listData.Error(iclus) = msg;
			continue
		end
	end

	% Get data for rest of table 
	[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, ...
		near_CF, oboe, bassoon, other] = finddata(params, CF);

	% If F0 isn't equal to 100 or 200 
	if any(bin | contra) && ~any(F0_100 | F0_200)
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

	listData.("Oboe")(iclus) = any(oboe);
	listData.("Bassoon")(iclus) = any(bassoon);
	listData.("Other")(iclus) = any(other);

end

%% Save spreadsheet

writetable(listData,'TimbreSessions_AllwithRVF.xlsx')


%% Functions 

function [bin, contra, F0_100, F0_200, level_43, level_63, level_73, ...
	level_83, even, near_CF, oboe, bassoon, other] = finddata(population, CF)

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

	% Find oboe 
	oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
		contains(p.list(1).wav_file, 'Oboe')), population);

	% Find bassoon
	bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
		contains(p.list(1).wav_file, 'Bassoon')), population);

	% Find other
	other = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), population);

end