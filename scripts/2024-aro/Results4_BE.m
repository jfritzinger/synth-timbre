%% BE_ModelPredictions
% J. Fritzinger, updated 9/25/23
%
% This script loads in all synthetic timbre data that has been saved and
% plots model responses (SFIE or broad inhibition model) to each session. 

clear all

%% Load in saved data

% Name of file
base = getPaths();
fpath = 'data/aro-2024';
filename = '2023-09-07_104117_SynthTimbre';

% Load Data
load(fullfile(base, fpath, [filename '.mat'])); % loads in 'population'

% Deletes any empty cells from the dataset
empty_cells = cellfun(@(p) isempty(p), population);
population(empty_cells) = [];
num_sessions = length(population);

%% Colors

data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};

%% Plot data and models


for isesh = 1:num_sessions

	% Get parameters and data for each session
	session = population{isesh};
	params = session(:,1);
	datas = session(:,2);
	rabbit = params{1}.animal;
	session = params{1}.session;
	tetrode = params{1}.tetrode;
	neuron = params{1}.neuron;

	% Find binaural/contra datasets
	is_contra = cellfun(@(p)isfield(p,'binmode')&&p.binmode ...
		== 1,params); % Sort data into binaural and contra
	is_ipsi = cellfun(@(p)isfield(p,'binmode')&&p.binmode ...
		== 0,params); % Sort data ipsi
	binmodes = {'Binaural', 'Contra', 'Ipsi'};
	binmode = 'Binaural';
	BIN = find(ismember(binmodes, binmode));
	params_bin = params(~is_contra);
	datas_bin = datas(~is_contra);

	has_mtf = cellfun(@(p)~isempty(p)&&strcmp(p.type,'typMTFN')&&...
		(p.all_mdepths(1) == 0), params_bin);
	has_rm = cellfun(@(p)~isempty(p)&&strcmp(p.type,'type=RM'), params_bin);
	has_sc = cellfun(@(p)~isempty(p)&&strcmp(p.type,'SPEC_slide')&&...
        strcmp(p.SPEC_slide_type,'Spectral_Centroid')&&...
        p.binmode==2&&mod(p.fpeak_mid, 200)==0&&p.Delta_F==200,params_bin);

	% Get MTF type
	if any(has_mtf)
		data_MTF = datas_bin{has_mtf};
		MTF_shape = data_MTF.MTF_shape;
		if contains(MTF_shape, 'H')
			MTF = 3;
			if strcmp(data_MTF.at_100, 'BE')
				MTF = 3;
				BMF = data_MTF.BMF;
			elseif strcmp(data_MTF.at_100, 'BS')
				MTF = 4;
			end
		elseif strcmp(MTF_shape, 'BE')
			MTF = 1;
			BMF = data_MTF.BMF;
		elseif strcmp(MTF_shape, 'BS')
			MTF = 2;
		else
			MTF = 5;
		end
	end

	% Only interested in BE units 
	if MTF ~= 1
		continue
	end

	if any(has_rm)
		data_RM = datas_bin{has_rm};
		CF = data_RM.CF;
	else
		continue
	end

    % Plot Spectral Centroid
    clear label
    data = datas_bin(has_sc);
	param_ST = params_bin(has_sc);
    num_DSIDs = length(data);
	if isempty(data) || isempty(data{1})
		continue
	end

	figure('Position', [460,427,1050,320])
    tiledlayout(1, 2,'TileSpacing','None')
    label_ind = 1;
    nexttile
    for ind = 1:num_DSIDs

        if param_ST{ind}.spl == 40 || param_ST{ind}.spl == 43
            color = data_colors{1};
        elseif param_ST{ind}.spl == 60 || param_ST{ind}.spl == 63
            color = data_colors{2};
        elseif param_ST{ind}.spl == 70 || param_ST{ind}.spl == 73
            color = data_colors{3};
        elseif param_ST{ind}.spl == 80 || param_ST{ind}.spl == 83
            color = data_colors{4};
        end
        DSID = param_ST{ind}.dsid;

        [fpeaks,~,fpeaksi] = unique([param_ST{ind}.list.fpeak].');
        num_fpeaks = length(fpeaks);
        dur = param_ST{ind}.dur/1000; % stimulus duration in seconds.
        rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];

		% Plot Data
		hold on
		h = errorbar(data{ind}.fpeaks,data{ind}.rate,data{ind}.rate_std/(sqrt(param_ST{ind}.nrep)),'LineWidth', 1.5, 'Color', color);
		if mod(param_ST{ind}.spl, 10)==0
			label(label_ind) = {[num2str(param_ST{ind}.spl+3) ' dB SPL']};
		else
			label(label_ind) = {[num2str(param_ST{ind}.spl) ' dB SPL']};
		end
		label_ind = label_ind+1;

    end

    % Plot Data
	%ylimits = [param_ST{ind}.fpeaks(1) param_ST{ind}.fpeaks(end)];
    xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    label(label_ind) = {'Estimated CF'};
    
    xlabel('Spectral Peak Frequency (Hz)')
    xlim([fpeaks(1) fpeaks(end)]);
    ylabel('Avg. Rate (sp/s)')
    %ylim(ylimits)
    set(gca,'FontSize',20)
	hLegend = legend(label, 'location', 'northeast', 'FontSize',14);
	hLegend.ItemTokenSize = [8,8];
    box on
    grid on
    title(['Single-Unit'], 'FontSize', 28);

    % Plot Model
    clear label
	clear stim_params AN SFIE
    nexttile
    label_ind = 1;
    for ind = 1:num_DSIDs

        % Parameters

		CS_param_names = {'Inhibitory strength low', 'Inhibitory strength high', 'Off-CF delay'}; % Parameters
		CS_param = [0.35, 0.35, 0.003];
		model_params.type = 'Lateral Model';
		model_params.config_type = 'BS inhibited by off-CF BS';
		model_params.lateral_CF = [CF/2 CF CF*2];
		model_params.CFs = model_params.lateral_CF;
		model_params.CF_range = model_params.CFs(2);

        stim_params.type = param_ST{ind}.type;
        stim_params.Fc = param_ST{ind}.fpeak_mid;
        stim_params.fpeak_mid = stim_params.Fc;
        stim_params.F0 = param_ST{ind}.Delta_F;
        stim_params.Fs = 100e3; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):
        stim_params.dur = param_ST{ind}.dur/1000;          % sec, stimulus duration
        stim_params.ramp_dur = param_ST{ind}.ramp_dur;     % sec, ramp duration
        stim_params.steps = param_ST{ind}.stp_otc;         % number of sliding stimuli
        stim_params.g = 24;
        stim_params.num_harms = param_ST{ind}.num_harms;
		stim_params.mnrep = 1;
		stim_params.physio = 0;
		stim_params.Delta_F = 200;
        if param_ST{ind}.spl < 100
            stim_params.spl = param_ST{ind}.spl;
        else
            stim_params.spl = 100;
        end

        model_params.species = 1;
        model_params.steps = stim_params.steps;
        %model_params.CF = CF;
        %model_params.BMF = BMF;
		model_params.CF_range = CF;
		model_params.num_CFs = 1;
		%model_params.CFs = CF;
		model_params.nAN_fibers_per_CF = 5;
		model_params.cohc = 1; % (0-1 where 1 is normal)
		model_params.cihc = 1; % (0-1 where 1 is normal)
		model_params.nrep = 3; % how many times to run the AN model
		model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
		model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
		model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
		model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
		% model_params.type = MTF_shape;
		model_params.BMF = 100;

        [stim_params] = generate_spectralcentroid(stim_params);
		model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
		%AN = modelAN(stim_params, model_params); % HSR for IC input
		%SFIE = wrapperIC(AN.an_sout, stim_params, model_params); % SFIE output

		AN = modelLateralAN(stim_params, model_params);
		SFIE = modelLateralSFIE(stim_params, model_params, AN.an_sout, AN.an_sout_lo, AN.an_sout_hi,'CS_params', CS_param);

        % Plot average IC response
		[avBE, ~] = plotSyntheticTimbre(stim_params, SFIE.avIC, 0);
        hold on
        plot(fpeaks,avBE, 'linewidth', 2, 'Color', model_colors{ind});
        
        % Calculate variance explained by the model
        R_int = corrcoef(data{ind}.rate,avBE).^2;
        R(isesh, ind) = R_int(1, 2);
        %[hat_r2er_BS(ind), r2_BS] = r2er_n2m2(avBS(ind, :)*1.7, BS.rate_matrix{ind});
        disp(['Variance explained = ' num2str(R(isesh, ind))]);

        label(label_ind) = {[num2str(stim_params.spl) ' dB SPL, R^2 = ' num2str(round(R(isesh, ind), 2))]};
        label_ind = label_ind+1;

    end

    % Plots
    xline(CF, '--', 'Color', [0.4 0.4 0.4],  'LineWidth', 2);
    label(label_ind) = {'Estimated CF'};
    
    xlim([fpeaks(1) fpeaks(end)]);
    %ylim(ylimits)
    yticklabels([])
    set(gca,'FontSize',20)
	hLegend = legend(label, 'location', 'northeast', 'FontSize',14);
	hLegend.ItemTokenSize = [8,8];
    xlabel('Spectral Peak Frequency (Hz)')
    title('Model Response', 'FontSize',28)
    box on
    grid on


end