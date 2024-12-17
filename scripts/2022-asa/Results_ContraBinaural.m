%% Results, Contra vs Binaural

% Single-unit BS that have contra/binaural
% R025S520_TT1N1 CF = 5000;
% R025S552_TT1N1 CF = 2000;
% R025S562_TT1N1 CF = 1500;
% R025S552_TT2N1 CF = 2000;
% R025S562_TT2N1 CF = 1500;
% R025S567_TT2N1 CF = 1300;
% R025S562_TT3N1 CF = 2300;
% R025S593_TT3N CF = 2000;
% R025S531_TT4N1 CF = 1200;
close all
clear all

%% Colors

RM_colors = {'#20116B', '#5E50A9', '#A49BD0'};
bin_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
contra_colors = {'#5BAEE8', '#1267B9', '#13357D', '#060650'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};

%% Load in datafiles

base = getPaths();
fpath = 'data/asa-2022';

%% Plot
% Peakier: [1 5 11]
% Other: [13 3 8]
for ii = [1, 3, 13, 8]

    switch ii
        case 1 % All levels, clear bin peakier than contra
            rabbit = 'R024';
            name = 'R024S478_TT2N1';
            CF = 1300;
            ylimits = [0 90];
        case 2 % All levels, no clear peaks
            rabbit = 'R025';
            name = 'R025S520_TT1N1';
            CF = 5000;
            ylimits = [0 35];
        case 3 % Yes, All levels, complicated but cool - mostly sharper except at 40dB
            rabbit = 'R025';
            name = 'R025S552_TT1N1';
            CF = 2000;
            ylimits = [0 100];
        case 4 % All levels, more complicated, not clear "peaks"
            rabbit = 'R025';
            name = 'R025S562_TT1N1';
            CF = 1500;
            ylimits = [0 110];
        case 5 % Yes, All levels, more complicated response changes
            rabbit = 'R025';
            name = 'R025S552_TT2N1';
            CF = 2000;
            ylimits = [0 110];
        case 6 % Yes, All levels, similar to case 8
            rabbit = 'R025';
            name = 'R025S562_TT2N1';
            CF = 1500;
            ylimits = [0 70];
        case 7 % Only have 60dB contra, seems similar peakiness
            rabbit = 'R025';
            name = 'R025S567_TT2N1';
            CF = 1300;
            ylimits = [0 70];
        case 8 % Yes, All levels, contra does not show peaks until 80dB, bin more sensitive
            rabbit = 'R025';
            name = 'R025S562_TT3N1';
            CF = 2300;
            ylimits = [0 80];
        case 9 % Contra is only 70dB, bin is sharper
            rabbit = 'R025';
            name = 'R025S593_TT3N1';
            CF = 2000;
            ylimits = [0 80];
        case 10 % two levels (70, 80) - bin is peakier than contra, but not best example
            rabbit = 'R025';
            name = 'R025S531_TT4N1';
            CF = 1200;
            ylimits = [0 50];
        case 11 % Yes, Has 3 levels, peaks are messy but it does show bin peakier than contra
            rabbit = 'R024';
            name = 'R024S460_TT2N1';
            CF = 1000;
            ylimits = [0 35];
            %BMF = 51; both
        case 12 % No, Only one contra, 60dB, the bin is sharper than contra
            rabbit = 'R024';
            name = 'R024S473_TT2N1';
            CF = 1150;
            ylimits = [0 120];
            %BMF = 122; both
        case 13 % Yes, Could be a good example after plotting other session too
            rabbit = 'R024';
            name = 'R024S484_TT2N1';
            name2 = 'R024S485_TT2N1';
            CF = 1500;
            ylimits = [0 90];
            %BMF = 40; both
    end

    % Figure Properties
    figure('position', [680,558,450, 800]); % left bottom width height
    tiledlayout(2,1,'TileSpacing','None');

    %% Plot binaural
    datafile = fullfile(base, fpath, rabbit, [name '.mat']);
    data_full = load(datafile, 'saved_data');
    data_full = data_full.saved_data;
    has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
        strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
        d.param.binmode == 2&&mod(d.param.fpeak_mid, 200)==0,data_full);

    % Plot Spectral Centroid
    clear label
    nexttile
    data = data_full(has_sc);
    num_DSIDs = length(data);

    label_ind = 1;
    for ind = 1:num_DSIDs
        DSID = data{ind}.param.dsid;
        if data{ind}.param.spl == 40
            color = bin_colors{1};
        elseif data{ind}.param.spl == 60
            color = bin_colors{2};
        elseif data{ind}.param.spl == 70
            color = bin_colors{3};
        elseif data{ind}.param.spl == 80
            color = bin_colors{4};
        end
        [fpeaks,~,fpeaksi] = unique([data{ind}.param.list.fpeak].');
        num_fpeaks = length(fpeaks);
        dur = data{ind}.param.dur/1000; % stimulus duration in seconds.

        rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
        spike_rates = data{ind}.spike_info.num_spikes_delayed/...
            (dur - data{ind}.spike_info.onsetWin/1000);

        if length(fpeaksi) == length(spike_rates)
            [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);

            % Plot
            hold on
            h = errorbar(fpeaks,rate,rate_std/(sqrt(data{ind}.param.nrep)),'LineWidth', 1.5, 'Color', color);
            label(label_ind) = {[num2str(data{ind}.param.spl+3) ' dB SPL']};
            label_ind = label_ind+1;
        end
    end
    xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    label(label_ind) = {'CF'};
    legend(label, 'location', 'best')
    %xlabel('Spectral Peak Frequency (Hz)')
    xticklabels([])
    xlim([fpeaks(1) fpeaks(end)]);
    ylabel('Avg. Rate (sp/s)')
    ylim(ylimits)
    set(gca,'FontSize',20)
    box on
    grid on
    %title([name ' Binaural']);

    %% Plot contra

    clear label
    has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
        strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
        d.param.binmode == 1,data_full);

    if ii == 13
        datacontra_full = load(fullfile(base, fpath, rabbit, [name2 '.mat']), 'saved_data');
        datacontra_full = datacontra_full.saved_data;
        has_sc_contra = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
            strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
            d.param.binmode == 1,datacontra_full);
        data_contra = datacontra_full(has_sc_contra);
    end

    nexttile
    data = data_full(has_sc);
    num_DSIDs = length(data);
    label_ind = 1;
    for ind = 1:num_DSIDs
        DSID = data{ind}.param.dsid;

        if data{ind}.param.spl == 40
            color = contra_colors{1};
        elseif data{ind}.param.spl == 60
            color = contra_colors{2};
        elseif data{ind}.param.spl == 70
            color = contra_colors{3};
        elseif data{ind}.param.spl == 80
            color = contra_colors{4};
        end
        [fpeaks,~,fpeaksi] = unique([data{ind}.param.list.fpeak].');
        num_fpeaks = length(fpeaks);
        dur = data{ind}.param.dur/1000; % stimulus duration in seconds.

        rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
        spike_rates = data{ind}.spike_info.num_spikes_delayed/...
            (dur - data{ind}.spike_info.onsetWin/1000);

        if length(fpeaksi) == length(spike_rates)
            [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);

            % Plot
            hold on
            h = errorbar(fpeaks,rate,rate_std/(sqrt(data{ind}.param.nrep)),'LineWidth', 1.5, 'Color', color);
            label(label_ind) = {[num2str(data{ind}.param.spl) ' dB SPL']};
            label_ind = label_ind+1;
        end
    end

    if ii == 13
        for ind = 1:length(data_contra)
            if data_contra{ind}.param.spl == 40
                color = contra_colors{1};
            elseif data_contra{ind}.param.spl == 60
                color = contra_colors{2};
            elseif data_contra{ind}.param.spl == 70
                color = contra_colors{3};
            elseif data_contra{ind}.param.spl == 80
                color = contra_colors{4};
            end
            DSID = data_contra{ind}.param.dsid;

            [fpeaks,~,fpeaksi] = unique([data_contra{ind}.param.list.fpeak].');
            num_fpeaks = length(fpeaks);
            dur = data_contra{ind}.param.dur/1000; % stimulus duration in seconds.

            rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
            spike_rates = data_contra{ind}.spike_info.num_spikes_delayed/...
                (dur - data_contra{ind}.spike_info.onsetWin/1000);

            if length(fpeaksi) == length(spike_rates)
                [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);

                % Plot
                hold on
                h = errorbar(fpeaks,rate,rate_std/(sqrt(data_contra{ind}.param.nrep)),'LineWidth', 1.5, 'Color', color);
                label(label_ind) = {[num2str(data_contra{ind}.param.spl) ' dB SPL']};
                label_ind = label_ind+1;
            end
        end
    end

    xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    label(label_ind) = {'CF'};
    legend(label, 'location', 'best')
    xlabel('Spectral Peak Frequency (Hz)')
    xlim([fpeaks(1) fpeaks(end)]);
    ylim(ylimits)
    ylabel('Avg. Rate (sp/s)')
    %yticklabels([])
    set(gca,'FontSize',20)
    box on
    grid on
    %title('Contralateral');

end
