%% Fig0_SavingData.m
% J. Fritzinger, updated 10/23/23
%
% Saves all data for WB-TIN paper
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

for isesh = 120:123

	% Name of file to save
	filename = sessions.Putative_Units{isesh};

	%% Load in each session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create struct for each dataset
	dataset_types = {'Binaural Response', 'Response Map', 'MTFN', 'SCHR', 'STRF',...
		'ST_43dB', 'ST_63dB', 'ST_73dB', 'ST_83dB', ...
		'ST_43dB_100', 'ST_63dB_100', 'ST_83dB_100', ...
		'Oboe', 'Bassoon', 'Other', 'RVF'};
	num_dataset_types = length(dataset_types);
	data = cell(num_dataset_types, 2);

	for ind = 10:34 

		% Get the path to each session
		full_name1 = sessions{isesh,ind};
		full_name = full_name1{1};

        if isempty(full_name)
            full_name1 = sessions{isesh,10};
    		full_name = full_name1{1};
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

            % Check if spreadsheet matches data
            [bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, ...
        		oboe, bassoon, other] = finddata(params);
            switch ind
                case 17
                    ST = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==2,params);
                case 18
                    ST = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==1,params);
                case 19
                    ST = cellfun(@(p)strcmp(p.type,'SCHR'),params);
                case 21
                    ST = any(bin & F0_200 & level_43 & even);
                case 22
    		        ST = any(bin & F0_200 & level_63 & even);
                case 23
    		        ST = any(bin & F0_200 & level_73 & even);
                case 24
    		        ST = any(bin & F0_200 & level_83 & even);
                case 25
        	        ST = any(bin & F0_100 & level_43 & even);
                case 26
                    ST = any(bin & F0_100 & level_63 & even);
                case 27
                    ST = any(bin & F0_100 & level_83 & even);
                case 28
    	            ST = any(contra & F0_200 & level_43 & even);
                case 29
    		        ST = any(contra & F0_200 & level_63 & even);
                case 30
                    ST = any(contra & F0_200 & level_73 & even);
                case 31
            		ST = any(contra & F0_200 & level_83 & even);
                case 32
            	    ST = any(oboe);
                case 33
            	    ST = any(bassoon);
                case 34
            	    ST = any(other);
            end
            if ind <17
                ST = 0;
            end
            if sum(ST)>0
               %error([full_name ': Forgot data for index ' num2str(ind)])
            end
            continue
        end

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

		%% Organize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Find sessions of interest
		putative = sessions.Putative_Units{isesh};
		for ids = 1:length(params)
			params{ids}.putative = putative;
			params{ids}.tetrode = tetrode;
			params{ids}.neuron = neuron;
		end

		% Check if this session has an spl shift
		spl_shift = checkIfSPLShift(rabbit, session);

		switch ind
			case 10 %%%%%% Binaural noise response
				has_bin = cellfun(@(p)strcmp(p.type,'char_spl'), params);
				data_binnoise = params{has_bin};
				data_binnoise.spls = [data_binnoise.spls(1:2)+spl_shift data_binnoise.spls(3)];
				for ilist = 1:length(data_binnoise.list)
					data_binnoise.list(ilist).spl = data_binnoise.list(ilist).spl+spl_shift;
				end
				data_binnoise.stim = get_stims(data_binnoise, stims);
				data_binnoise.cluster = cluster;
				data{1,2} = data_binnoise;

			case 11 % ITD
				% Not saving
			case 12 % ILD
				% Not saving
			case 13 %%%%%% RM binaural
				has_rm = cellfun(@(p)strcmp(p.type,'type=RM')&&p.binmode==2, params);
				if any(has_rm)
					data_rm = params{has_rm};
					data_rm.spls = [data_rm.spls(1:2)+spl_shift data_rm.spls(3)];
					data_rm.all_spls = data_rm.all_spls+spl_shift;
					for ilist = 1:length(data_rm.list)
						data_rm.list(ilist).spl = data_rm.list(ilist).spl+spl_shift;
					end
					data_rm.stims = get_stims(data_rm, stims);
					data_rm.cluster = cluster;
					data{2,2} = data_rm;
				end

			case 14 %%%%%% RM contra
				has_rm = cellfun(@(p)strcmp(p.type,'type=RM')&&p.binmode==1, params);
				if any(has_rm)
					data_rm = params{has_rm};
					data_rm.spls = [data_rm.spls(1:2)+spl_shift data_rm.spls(3)];
					data_rm.all_spls = data_rm.all_spls+spl_shift;
					for ilist = 1:length(data_rm.list)
						data_rm.list(ilist).spl = data_rm.list(ilist).spl+spl_shift;
					end
					data_rm.stims = get_stims(data_rm, stims);
					data_rm.cluster = cluster;
					data{2,1} = data_rm;
				end

			case 15 %%%%%% MTFN binaural
				has_mtf = cellfun(@(p)strcmp(p.type,'typMTFN')&&(p.all_mdepths(1) == 0)&&...
					p.binmode==2, params);
				if any(has_mtf)
					data_mtf = params{has_mtf};
					data_mtf.spl = data_mtf.spl+spl_shift;
					data_mtf.stims = get_stims(data_mtf, stims);
					data_mtf.cluster = cluster;
					data{3,2} = data_mtf;
				end

			case 16 %%%%%% MTFN contra
				has_mtf = cellfun(@(p)strcmp(p.type,'typMTFN')&&(p.all_mdepths(1) == 0)&&...
					p.binmode==1, params);
				if any(has_mtf)
					data_mtf = params{has_mtf};
					data_mtf.spl = data_mtf.spl+spl_shift;
					data_mtf.stims = get_stims(data_mtf, stims);
					data_mtf.cluster = cluster;
					data{3,1} = data_mtf;
				end

			case 17 %%%%%% STRF binaural
				has_strf = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==2,params);
				if any(has_strf)
					data_strf = params{has_strf};
					data_strf.spl = data_strf.spl + spl_shift;
					data_strf.stims = get_stims(data_strf, stims);
					data_strf.cluster = cluster;
					data{4,2} = data_strf;
				end

			case 18 %%%%%% STRF contra
				has_strf = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==1,params);
				if any(has_strf)
					data_strf = params{has_strf};
					data_strf.spl = data_strf.spl + spl_shift;
					data_strf.stims = get_stims(data_strf, stims);
					data_strf.cluster = cluster;
					data{4,1} = data_strf;
				end

			case 19 %%%%%% SCHR
				has_schr = cellfun(@(p)strcmp(p.type,'SCHR'),params);
				if any(has_schr)
					data_schr = params{has_schr};
					data_schr.stimdB = data_schr.stimdB+spl_shift;
					data_schr.stim = get_stims(data_schr, stims);
					data_schr.cluster = cluster;
					data{5,2} = data_schr;
				end
			case 20 
				has_rvf = cellfun(@(p)strcmp(p.type,'RVF'),params);
				if any(has_rvf)
					data_rvf = params{has_rvf};
					data_rvf.stim = get_stims(data_rvf, stims);
					data_rvf.cluster = cluster;
					data{5,1} = data_rvf;
				end
			case 21 %%%%%% Synth Timbre, 43 dB SPL, binaural (row 6)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==40||p.spl==43)&&p.binmode==2&&p.Delta_F==200&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{6,2} = data_sc;
				end

			case 22 %%%%%% Synth Timbre, 63 dB SPL, binaural (row 7)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==60||p.spl==63)&&p.binmode==2&&p.Delta_F==200&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{7,2} = data_sc;
				end

			case 23 %%%%%% Synth Timbre, 73 dB SPL, binaural (row 8)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==70||p.spl==73)&&p.binmode==2&&p.Delta_F==200&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{8,2} = data_sc;
				end

			case 24 %%%%%% Synth Timbre, 83 dB SPL, binaural (row 9)
				has_sc_83 = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==80||p.spl==83)&&p.binmode==2&&p.Delta_F==200&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc_83)
				if any(has_sc_83)
					data_sc_83 = params{has_sc_83};
					data_sc_83.spl = data_sc_83.spl+spl_shift;
					data_sc_83.stims = get_stims(data_sc_83, stims);
					data_sc_83.cluster = cluster;
					data{9,2} = data_sc_83;
				end

			case 25 %%%%%% Synth Timbre, 43 dB SPL, 100Hz, bin (row 10)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==40||p.spl==43)&&p.binmode==2&&p.Delta_F==100&& mod(p.fpeak_mid, 100)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{10,2} = data_sc;
				end

			case 26 %%%%%% Synth Timbre, 63 dB SPL, 100Hz, bin (row 11)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==60||p.spl==63)&&p.binmode==2&&p.Delta_F==100&& mod(p.fpeak_mid, 100)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{11,2} = data_sc;
				end

			case 27 %%%%%% Synth Timbre, 83 dB SPL, 100Hz, bin (row 12)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==80||p.spl==83)&&p.binmode==2&&p.Delta_F==100&& mod(p.fpeak_mid, 100)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{12,2} = data_sc;
				end
			case 28 %%%%%% Synth Timbre, 43 dB SPL, contra (row 6)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==40||p.spl==43)&&p.binmode==1&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{6,1} = data_sc;
				end

			case 29 %%%%%% Synth Timbre, 63 dB SPL, contra (row 7)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==60||p.spl==63)&&p.binmode==1&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{7,1} = data_sc;
				end

			case 30 %%%%%% Synth Timbre, 73 dB SPL, contra (row 7)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==70||p.spl==73)&&p.binmode==1&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{8,1} = data_sc;
				end

			case 31 %%%%%% Synth Timbre, 83 dB SPL, contra (row 7)
				has_sc = cellfun(@(p)strcmp(p.type,'SPEC_slide')&&strcmp(p.SPEC_slide_type,'Spectral_Centroid')...
					&&(p.spl==80||p.spl==83)&&p.binmode==1&& mod(p.fpeak_mid, 200)==0,params);
				checkIfMultiple(has_sc)
				if any(has_sc)
					data_sc = params{has_sc};
					data_sc.spl = data_sc.spl+spl_shift;
					data_sc.stims = get_stims(data_sc, stims);
					data_sc.cluster = cluster;
					data{9,1} = data_sc;
				end
			case 32 %%%%%% Oboe
				has_oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
					(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
					contains(p.list(1).wav_file, 'Oboe')), params);
				if any(has_oboe)
					data_oboe = params{has_oboe};
					data_oboe.spl = data_oboe.signal_spls+spl_shift;
					data_oboe.stims = get_stims(data_oboe, stims);
					data_oboe.cluster = cluster;
					data{13,2} = data_oboe;
				end

			case 33 %%%%%% Bassoon
				has_bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
					(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
					contains(p.list(1).wav_file, 'Bassoon')), params);
				if any(has_bassoon)
					data_bassoon = params{has_bassoon};
					data_bassoon.spl = data_bassoon.signal_spls+spl_shift;
					data_bassoon.stims = get_stims(data_bassoon, stims);
					data_bassoon.cluster = cluster;
					data{14,2} = data_bassoon;
				end

			case 34 %%%%%% Other natural timbre
				has_nt = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), params);
				index = find(has_nt);

				if any(has_nt)
					for ii = 1:length(index)
						jj= index(ii);
						idata = 14+ii;

						data_nt = params{jj};
						data_nt.spl = data_nt.signal_spls+spl_shift;
						data_nt.stims = get_stims(data_nt, stims);
						data_nt.cluster = cluster;
						data{idata,2} = data_nt;
					end
				end
		end
	end

	% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	save(fullfile(path, filename), 'data')

end


%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to get the stimulus properties for each DSID
function [stim] = get_stims(data, stims)
ds = data.dsid;
stim = stims(ds);
end

% Function that checks to see if the session had the +3dB error in
% recording, and if so sends back the +3 dB value.
function spl_shift = splShift(rabbit, session)
spl_shift = 0;
switch rabbit
	case 24
		spl_shift = 3;
	case 25
		if session < 676
			spl_shift = 3;
		end
	case 27
		if session < 81
			spl_shift = 3;
		end
end
end


function checkIfMultiple(array)
	if sum(array)>1
		error('Two possible conditions to save, check!')
	elseif sum(array)==0
		error('Spreadsheet and data do not match!')
	end
end

function [bin, contra, F0_100, F0_200, level_43, level_63, level_73, ...
	level_83, even, oboe, bassoon, other] = finddata(population)

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
	%datasets = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
	%	strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);
	%near_CF = cellfun(@(p) min(abs(p.fpeak_mid-CF)), population(datasets));

	% Find oboe 
	oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
		contains(p.list(1).wav_file, 'Oboe')), population);

	% Find bassoon
	bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
		contains(p.list(1).wav_file, 'Bassoon')), population);

	% Find other
	other = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), population)...
        & ~bassoon & ~oboe;

end