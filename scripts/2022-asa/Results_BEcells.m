%% BE data

% Single-unit BS cells that are matched well to CF = Fpeak
% R025S467
% R025S520
% R025S529
% R025S388
% R025S474
%
% R24
% R024S418_TT2N1 CF = 4000; y = [0 40]; BMF = 354;
% R024S444_TT2N1 CF = 5000; y = [0 70]; BMF = 300;
%
% R27:
close all

%% Colors

data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};

%% Find filepaths

base = getPaths();
fpath = 'data/asa-2022';

%% Plot data and models

% 3, 8
% 4, 5, 2
% Figure Properties
BMF = 100;

for ii = [3, 8, 4, 5, 2]
    switch ii
        case 1 % Shows decreasing response, doesn't match model (model has dips)
            rabbit = 'R025';
            session = 'R025S474_TT1N2';
            CF = 4600;
            ylimits = [0 40];
            %BMF = 82;
        case 2 % Shows peaks in responses, model response doesn't show much (40, 60 dB)
            rabbit = 'R025';
            session = 'R025S529_TT3N1';
            CF = 5000;
            ylimits = [0 90];
            %BMF = 200;
        case 3 % Slightly matched in model, need to use only 1000 (not 1100), one level 70dB
            rabbit = 'R025';
            session = 'R025S388_TT4N2';
            CF = 950;
            ylimits = [0 140];
            %BMF = 66;
        case 4 % All levels, shows decreasing response not matched in model 
            rabbit = 'R025';
            session = 'R025S520_TT2N1';
            CF = 4600;
            ylimits = [0 70];
            %BMF = 65;
        case 5 % Off center by a lot, show flat response , all levels 
            rabbit = 'R025';
            session = 'R025S467_TT2N1';
            CF = 3000;
            ylimits = [0 70];
            %BMF = 86;
        case 6 % All 3 levels, 40dB slightly matches but otherwise flat response and model shows dips 
            rabbit = 'R024';
            session = 'R024S418_TT2N1';
            CF = 4000;
            ylimits = [0 40];
            %BMF = 354;
        case 7 % only one level, nothing remarkable to say about this one - mostly flat 
            rabbit = 'R024';
            session = 'R024S444_TT2N1';
            CF = 5000;
            ylimits = [0 70];
            %BMF = 300;
        case 8 % Clearly shows peaks when dips are expected, all levels 
            rabbit = 'R027';
            session = 'R027S037_TT3N1';
            CF = 1740;
            ylimits = [0 70];
            %BMF = 62;
    end
    disp(session)

    figure('Position', [460,427,1050,320])

    % Load data
    datafile = fullfile(base, fpath, rabbit, session);
    data_full = load(datafile, 'saved_data');
    data_full = data_full.saved_data;
    has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
        strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
        d.param.binmode==2&&mod(d.param.fpeak_mid, 200)==0,data_full);

    % Plot Spectral Centroid
    clear label
    data = data_full(has_sc);
    num_DSIDs = length(data);
    label_ind = 1;
    nexttile
    for ind = 1:num_DSIDs
        if data{ind}.param.spl == 40
            color = data_colors{1};
        elseif data{ind}.param.spl == 60
            color = data_colors{2};
        elseif data{ind}.param.spl == 70
            color = data_colors{3};
        elseif data{ind}.param.spl == 80
            color = data_colors{4};
        end
        DSID = data{ind}.param.dsid;

        [fpeaks,~,fpeaksi] = unique([data{ind}.param.list.fpeak].');
        num_fpeaks = length(fpeaks);
        dur = data{ind}.param.dur/1000; % stimulus duration in seconds.

        rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
        spike_rates = data{ind}.spike_info.num_spikes_delayed/...
            (dur - data{ind}.spike_info.onsetWin/1000);

        if length(fpeaksi) == length(spike_rates)
            [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);

            % Plot Data
            hold on
            h = errorbar(fpeaks,rate,rate_std/(sqrt(data{ind}.param.nrep)),'LineWidth', 1.5, 'Color', color);
            label(label_ind) = {[num2str(data{ind}.param.spl) ' dB SPL']};
            label_ind = label_ind+1;
        end
    end

    % Plot Data
    xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    label(label_ind) = {'Estimated CF'};
    legend(label, 'location', 'northeast')
    xlabel('Spectral Peak Frequency (Hz)')
    xlim([fpeaks(1) fpeaks(end)]);
    ylabel('Avg. Rate (sp/s)')
    ylim(ylimits)
    set(gca,'FontSize',20)
    box on
    grid on
    if ii ~= 8 && ii ~= 5 && ii ~= 2
        title('Single-Unit', 'FontSize', 28);
    end

    % Plot Model
    clear label
    nexttile
    label_ind = 1;
    for ind = 1:num_DSIDs

        % Parameters
        params.type = data{ind}.param.type;
        params.Fc = data{ind}.param.fpeak_mid;
        params.fpeak_mid = params.Fc;
        params.F0 = data{ind}.param.Delta_F;
        params.Fs = 100e3; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):
        params.dur = data{ind}.param.dur/1000;          % sec, stimulus duration
        params.ramp_dur = data{ind}.param.ramp_dur;     % sec, ramp duration
        params.steps = data{ind}.param.stp_otc;         % number of sliding stimuli
        params.G = 24;
        params.num_harms = data{ind}.param.num_harms;
		params.freq_lo = params.fpeak_mid - (params.num_harms-1)/2 * params.Delta_F; % three components below (try 2.5 7/30/20 LHC)
		params.freq_hi = params.fpeak_mid + (params.num_harms-1)/2 * params.Delta_F; % three components above ( ditto )
		if params.freq_lo < 200 % limited to 200 Hz on upper end
			params.freq_lo = 200;
		end
		if params.freq_hi > 20000 % limited to 20 kHz on upper end
			params.freq_hi = 20000;
		end

        if data{ind}.param.spl < 100
            params.stimdB = data{ind}.param.spl;
        else
            params.stimdB = 100;
        end

        model_params.species = 1;
        model_params.steps = params.steps;
        model_params.CF = CF;
        model_params.BMF = BMF;
        [stim, model_params.stim_list, model_params.fpeaks] = generate_ST(params);
        [~, ~, avBE, stdBE, avBS, stdBS, ~] = modelSingleCell(stim, model_params, params);

        % Plot average IC response
        hold on
        plot(fpeaks,avBE, 'linewidth', 2, 'Color', model_colors{ind});
        
        % Calculate variance explained by the model
        R_int = corrcoef(rate,avBE).^2;
        R(ii, ind) = R_int(1, 2);
        %[hat_r2er_BS(ind), r2_BS] = r2er_n2m2(avBS(ind, :)*1.7, BS.rate_matrix{ind});
        disp(['Variance explained = ' num2str(R(ii, ind))]);

        label(label_ind) = {[num2str(params.stimdB) ' dB SPL, R^2 = ' num2str(round(R(ii, ind), 2))]};
        label_ind = label_ind+1;

    end

    % Plots
    xline(CF, '--', 'Color', [0.4 0.4 0.4],  'LineWidth', 2);
    label(label_ind) = {'Estimated CF'};
    legend(label, 'location', 'northeast')
    xlim([fpeaks(1) fpeaks(end)]);
    ylim(ylimits)
    yticklabels([])
    set(gca,'FontSize',20)
    xlabel('Spectral Peak Frequency (Hz)')
    if ii ~= 8  && ii ~= 5 && ii ~= 2
        title('Model Response', 'FontSize',28)
    end
    box on
    grid on


end


