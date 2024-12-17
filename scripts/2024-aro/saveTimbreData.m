%% Save timbre data for ARO
% J. Fritzinger
clear
%% Load in spreadsheet 

spreadsheet_name = 'SynthSessions.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

base = getPaths();
fpath = 'data/aro-2024';

[userid, base_dir, ~, report_path, data_path] = findPaths();

%% Load and save data 

for iclus = [20, 80, 88, 118]

	% Get the path to each session
	session = sessions.Session{iclus};
	rab_num = num2str(sessions.Rabbit(iclus));
	if ismac
		session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
	else							
		session_dir_name = base_dir{contains(base_dir, rab_num)};
	end
	session_dir = fullfile(session_dir_name, session);

	% Run post process
	[clusters, params, stims] = loadPhysiologySession(session_dir, session, userid);
	CF = sessions.CF(iclus);

	% Get data for tetrode/neuron of interest
	cluster = clusters([clusters.tetrode] == sessions.TT(iclus)); % Select tetrode
	cluster = cluster([cluster.neuron] == sessions.N(iclus)); % Select neuron

	% Find datasets
	has_bin = cellfun(@(p)strcmp(p.type,'char_spl'), params);
	has_charITD = cellfun(@(p)strcmp(p.type,'char_ITD'), params);
	has_charILD = cellfun(@(p)strcmp(p.type,'char_ILD'), params);
	has_rm = cellfun(@(p)strcmp(p.type,'type=RM'), params);
	has_mtf = cellfun(@(p)strcmp(p.type,'typMTFN')&&(p.all_mdepths(1) == 0), params);
	has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid'),params);
	has_strf = cellfun(@(p) strcmp(p.type,'STRF'),params);
	has_nt = cellfun(@(p)strcmp(p.type,'Natural_Timbre'),params);

	% Set up data
	num_dsids = size(params, 1);
	data = cell(num_dsids, 1);

	% Binaural Response
	if any(has_bin)
		[fig, data{has_bin}] = plotPhysBIN(cluster, params(has_bin),stims);
		close
	end

	% Characterizing ITD
	if any(has_charITD)
		[fig, data{has_charITD}] = plotPhysCharITD(cluster, params(has_charITD),stims);
		close
	end

	% Characterizing ILD
	if any(has_charILD)
		[fig, data{has_charILD}] = plotPhysCharILD(cluster, params(has_charILD),stims);
		close
	end

	% RM
	if any(has_rm)
		num_rm = sum(has_rm);
		rm_ds = find(has_rm);
		for j = 1:num_rm
			[fig, data{rm_ds(j)}] = plotPhysRM(cluster, params(rm_ds(j)),stims, CF);
			close
		end
	end

	% MTF
	if any(has_mtf)
		num_mtf = sum(has_mtf);
		mtf_ds = find(has_mtf);
		for j = 1:num_mtf
			[~,~,~, fig, data{mtf_ds(j)}] = plotPhysMTFTTest(cluster, params(mtf_ds(j)),stims);
			close
		end
	end

	% STRF
	if any(has_strf)
		num_strf = sum(has_strf);
		strf_ds = find(has_strf);
		for j = 1:num_strf
			[params(strf_ds(j)), fig, data{strf_ds(j)}] = plotPhysSTRF(cluster,params(strf_ds(j)), stims);
			close
		end
	end

	% Avg Rate Synthetic Timbre
	% Saved average rate data
	if any(has_sc)
		[params(has_sc), ~, data(has_sc)] = plotPhysST(cluster, params(has_sc), stims, CF, data(has_sc), []);
		close
		[params(has_sc), ~, data(has_sc)] = plotPhysST_Rasters(cluster, params(has_sc), stims, 'PSTH', '', data(has_sc));
		close 
	end

	% Natural Timbre
	if any(has_nt)
		[params(has_nt), ~, data(has_nt)] = plotPhysNT(cluster, params(has_nt),stims, 0, data(has_nt));
		close 
	end


	filename = sprintf('%s_TT%d_N%d.mat', session, sessions.TT(iclus), sessions.N(iclus));

	%% Save data 
	save(fullfile(base, fpath, filename), 'params', 'data', 'cluster', 'stims')

end