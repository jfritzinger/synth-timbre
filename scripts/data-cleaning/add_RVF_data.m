%% add_RVF_data.m
clear

%% Load in PDF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

if ismac
	path = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/data/manuscript/neural_data';
else
	path = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre\data\manuscript\neural_data';
end

for isesh = 303 % 218, 256, 303 error

	filename = sessions.Putative_Units{isesh};
	full_name = sessions.RVF{isesh};
	rabbit = str2num(full_name(3:4));
	session = str2num(full_name(6:8));
	tetrode = str2num(full_name(12));
	neuron = str2num(full_name(15));
	[userid, base_dir, ~, report_path, data_path] = findPaths();
	session_name = sprintf('R%03dS%03d', rabbit, session);
	rab_num = num2str(rabbit);
	if ismac
		session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
	else
		session_dir_name = base_dir{contains(base_dir, rab_num)};
	end
	session_dir = fullfile(session_dir_name, session_name);

	% Post_process
	[clusters, params, stims] = loadPhysiologySession(session_dir, session_name, userid);

	% Get data for tetrode/neuron of interest
	cluster = clusters([clusters.tetrode] == tetrode); % Select tetrode
	cluster = cluster([cluster.neuron] == neuron); % Select neuron

	% Load in other
	load(fullfile(path, filename), 'data')

	% Add RVF
	has_rvf = cellfun(@(p)strcmp(p.type,'RVF'),params);
	if any(has_rvf)
		data_rvf = params{has_rvf};
		data_rvf.stim = get_stims(data_rvf, stims);
		data_rvf.cluster = cluster;
		data{5,1} = data_rvf;
	end

	% Resave
	save(fullfile(path, filename), 'data')
end

%% Function to get the stimulus properties for each DSID
function [stim] = get_stims(data, stims)
ds = data.dsid;
stim = stims(ds);
end