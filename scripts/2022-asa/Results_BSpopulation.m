%% Figure 8: Summary of many datasets
% J. Fritzinger, updated 2/23/2021
clear all
close all

% Single-unit BS cells that are matched well to CF = Fpeak
%
% R025:
% R025S438_TT1N1 CF = 2000; y = []; BMF = ;
% R025S442_TT1N1 CF = 2600;
% R025S474_TT1N1 CF = 4600;
% R025S520_TT1N1 CF = 5000;
% R025S552_TT1N1 CF = 2000;
% R025S562_TT1N1 CF = 1500;
% R025S372_TT2N2 CF = 1300;
% R025S552_TT2N1 CF = 2000;
% R025S562_TT2N1 CF = 1500;
% R025S567_TT2N1 CF = 1300;
% R025S442_TT3N1 CF = 2300;
% R025S562_TT3N1 CF = 2300;
% R025S593_TT3N CF = 2000;
% R025S531_TT4N1 CF = 1200;
%
% R24: No TT4
% R024S410_TT2N1 CF = 1740; y = [0 25]; BMF = 80; no contra
% R024S450_TT2N1 CF = 1760; y = [0 100]; BMF = 43; no contra
% R024S460_TT2N1 CF = 1000; y = [0 35]; BMF = 51; both
% R024S470_TT2N1 CF = 1170; y = [0 80]; BMF = 181; no contra
% R024S473_TT2N1 CF = 1150; y = [0 120]; BMF = 122; both
% R024S478_TT2N1 CF = 1300; y = [0 90]; BMF = 80; both
% R024S484_TT2N1 CF = 1500; y = [0 90]; BMF = 40; both also: R024S485, contra
% R024S402_TT1N2 CF = 400; y = [0 120]; BMF = 27; no contra
%
% R27:
% 
% 

% Hybrid:
% R024S407_TT2N1, shows peaks 
% R025S474_TT1N2
% R025S484_TT1N1
% R025S470_TT2N1
% R025S552_TT3N1
% R025S552_TT4N1
% R025S562_TT4N1


%% Load in datafiles

base = getPaths();
fpath = 'data/asa-2022';

pop_colors = {'#E4C9AE', '#E8B88A','#E8A361','#E8882C', '#C06A17', '#E28F3F', '#DE7511', ...
    '#A75B13', '#7E430B', '#3E2003', '#271402'};

%% Datasets
% 
% j = 1;
% spl = [40 60 80];
% for spl_ind = 1:13
%     data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C', '#000000'};
%     fontname = 'Arial';
%     set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
%     figure;
%     set(gcf, 'Position', [560,527,250,250])
%     set(gcf,'color','w');
% 
%     label = {};
%     label_ind = 1;
%     j = 1;
%     CFs = [];
%     rate_diffs = [];
%     for ii = 1:23
%         switch ii
%             case 1 % R25
%                 rabbit = 'R025';
%                 session = 'R025S438_TT1N1'; CF = 2000;
%             case 2
%                 rabbit = 'R025';
%                 session = 'R025S442_TT1N1'; CF = 2600;
%             case 3
%                 rabbit = 'R025';
%                 session = 'R025S474_TT1N1'; CF = 4600;
%             case 4
%                 rabbit = 'R025';
%                 session = 'R025S520_TT1N1'; CF = 5000;
%             case 5
%                 rabbit = 'R025';
%                 session = 'R025S552_TT1N1'; CF = 2000;
%             case 6
%                 rabbit = 'R025';
%                 session = 'R025S562_TT1N1'; CF = 1500;
%             case 7
%                 rabbit = 'R025';
%                 session = 'R025S372_TT2N2'; CF = 1300;
%             case 8 
%                 rabbit = 'R025';
%                 session = 'R025S552_TT2N1'; CF = 2000;
%             case 9
%                 rabbit = 'R025';
%                 session = 'R025S562_TT2N1'; CF = 1500;
%             case 10
%                 rabbit = 'R025';
%                 session = 'R025S567_TT2N1'; CF = 1300;
%             case 11
%                 rabbit = 'R025';
%                 session = 'R025S442_TT3N1'; CF = 2300;
%             case 12
%                 rabbit = 'R025';
%                 session = 'R025S562_TT3N1'; CF = 2300;
%             case 13
%                 rabbit = 'R025';
%                 session = 'R025S593_TT3N1'; CF = 2000;
%             case 14
%                 rabbit = 'R025';
%                 session = 'R025S531_TT4N1'; CF = 1200;
%             case 15
%                 rabbit = 'R024';
%                 session = 'R024S410_TT2N1'; CF = 1740; %y = [0 25]; BMF = 80; no contra
%             case 16
%                 rabbit = 'R024';
%                 session = 'R024S450_TT2N1'; CF = 1760; %y = [0 100]; BMF = 43; no contra
%             case 17
%                 rabbit = 'R024';
%                 session = 'R024S460_TT2N1'; CF = 1000; %y = [0 35]; BMF = 51; both
%             case 18
%                 rabbit = 'R024';
%                 session = 'R024S470_TT2N1'; CF = 1170; %y = [0 80]; BMF = 181; no contra
%             case 19
%                 rabbit = 'R024';
%                 session = 'R024S473_TT2N1'; CF = 1150; %y = [0 120]; BMF = 122; both
%             case 20
%                 rabbit = 'R024';
%                 session = 'R024S478_TT2N1'; CF = 1300; %y = [0 90]; BMF = 80; both
%             case 21
%                 rabbit = 'R024';
%                 session = 'R024S484_TT2N1'; CF = 1500; %y = [0 90]; BMF = 40; both also: R024S485, contra
%             case 22
%                 rabbit = 'R024';
%                 session = 'R024S402_TT1N2'; CF = 400; %y = [0 120]; BMF = 27; no contra
%             case 23
%                 rabbit = 'R027';
%                 session = 'R027S019_TT3N1'; CF = 1000;
%                 %session = 'R027S019_TT3N1'; CF = 1000;
%         end
% 
%         datafile = [fpath rabbit slash session '.mat'];
%         data_full = load(datafile, 'saved_data');
%         data_full = data_full.saved_data;
%         has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
%             strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
%             d.param.binmode == 2&&d.param.spl==spl(spl_ind),data_full);
% 
%         % Plot Spectral Centroid
%         data = data_full(has_sc);
%         num_DSIDs = length(data);
%         for ind = 1:num_DSIDs
%             DSID = data{ind}.param.dsid;
% 
%             [fpeaks,~,fpeaksi] = unique([data{ind}.param.list.fpeak].');
%             num_fpeaks = length(fpeaks);
%             dur = data{ind}.param.dur/1000; % stimulus duration in seconds.
% 
%             rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
%             spike_rates = data{ind}.spike_info.num_spikes_delayed/...
%                 (dur - data{ind}.spike_info.onsetWin/1000);
% 
%             if length(fpeaksi) == length(spike_rates)
%                 [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);
%                 CFs(j) = CF;
%                 
%                 % Run model
%                 % Parameters
%                 params.type = data{ind}.param.type;
%                 params.Fc = data{ind}.param.fpeak_mid;
%                 params.fpeak_mid = params.Fc;
%                 params.F0 = data{ind}.param.Delta_F;
%                 params.Fs = 100e3; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):
%                 params.dur = data{ind}.param.dur/1000;          % sec, stimulus duration
%                 params.ramp_dur = data{ind}.param.ramp_dur;     % sec, ramp duration
%                 params.steps = data{ind}.param.stp_otc;         % number of sliding stimuli
%                 params.G = 24;
%                 params.num_harms = data{ind}.param.num_harms;
%                 [params.freq_lo, params.freq_hi] = get_freq_limits(params.Fc, params.num_harms, params.F0);
%                 if data{ind}.param.spl < 100
%                     params.stimdB = data{ind}.param.spl;
%                 else
%                     params.stimdB = 100;
%                 end
% 
%                 model_params.species = 1; % 1 = cat 
%                 model_params.steps = params.steps;
%                 model_params.CF = CF;
%                 model_params.BMF = 100;
%                 [stim, model_params.stim_list, model_params.fpeaks] = generate_spectralcentroid(params);
%                 [~, ~, avBE, stdBE, avBS, stdBS, ~] = modelSingleCell(stim, model_params, params);
% 
%                 % Calculate variance explained by the model
%                 R_int = corrcoef(rate,avBS).^2;
%                 R(j) = R_int(1, 2);
%                 %[hat_r2er_BS(ind), r2_BS] = r2er_n2m2(avBS(ind, :)*1.7, BS.rate_matrix{ind});
%                 disp(['Variance explained = ' num2str(R(j))]);
% 
%                 j = j+1;
%             end
%         end  
%     end
%     scatter(CFs, R, 'filled')
%     r2{ind}.CFs = CFs;
%     r2{ind}.R = R;
%     xlabel('CF')
%     xlim([100 5000]);
%     ylabel('Peak Height (sp/s)')
%     ylim([0 1])
%     set(gca,'FontSize',16)
%     %set(gca, 'XScale', 'log')
%     box on
%     grid on
%     title([num2str(spl_ind) 'dB SPL']);
% end

%% Align to CF

fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
figure;
set(gcf, 'Position', [560,527,1100,500])
set(gcf,'color','w');
spl = [43, 63, 83];
spl2 = [40, 60, 80];

for spl_ind = 1:3
    label = {};
    j = 1;
    clear index CFs
    for ii = 1:23
        switch ii
            case 1 % R25
                rabbit = 'R025';
                session = 'R025S438_TT1N1'; CF = 2000;
            case 2
                rabbit = 'R025';
                session = 'R025S442_TT1N1'; CF = 2600;
            case 3
                rabbit = 'R025';
                session = 'R025S474_TT1N1'; CF = 4600;
            case 4
                rabbit = 'R025';
                session = 'R025S520_TT1N1'; CF = 5000;
            case 5
                rabbit = 'R025';
                session = 'R025S552_TT1N1'; CF = 2000;
            case 6
                rabbit = 'R025';
                session = 'R025S562_TT1N1'; CF = 1500;
            case 7
                rabbit = 'R025';
                session = 'R025S372_TT2N2'; CF = 1300;
            case 8 
                rabbit = 'R025';
                session = 'R025S552_TT2N1'; CF = 2000;
            case 9
                rabbit = 'R025';
                session = 'R025S562_TT2N1'; CF = 1500;
            case 10
                rabbit = 'R025';
                session = 'R025S567_TT2N1'; CF = 1300;
            case 11
                rabbit = 'R025';
                session = 'R025S442_TT3N1'; CF = 2300;
            case 12
                rabbit = 'R025';
                session = 'R025S562_TT3N1'; CF = 2300;
            case 13
                rabbit = 'R025';
                session = 'R025S593_TT3N1'; CF = 2000;
            case 14
                rabbit = 'R025';
                session = 'R025S531_TT4N1'; CF = 1200;
            case 15
                rabbit = 'R024';
                session = 'R024S410_TT2N1'; CF = 1740; %y = [0 25]; BMF = 80; no contra
            case 16
                rabbit = 'R024';
                session = 'R024S450_TT2N1'; CF = 1760; %y = [0 100]; BMF = 43; no contra
            case 17
                rabbit = 'R024';
                session = 'R024S460_TT2N1'; CF = 1000; %y = [0 35]; BMF = 51; both
            case 18
                rabbit = 'R024';
                session = 'R024S470_TT2N1'; CF = 1170; %y = [0 80]; BMF = 181; no contra
            case 19
                rabbit = 'R024';
                session = 'R024S473_TT2N1'; CF = 1150; %y = [0 120]; BMF = 122; both
            case 20
                rabbit = 'R024';
                session = 'R024S478_TT2N1'; CF = 1300; %y = [0 90]; BMF = 80; both
            case 21
                rabbit = 'R024';
                session = 'R024S484_TT2N1'; CF = 1500; %y = [0 90]; BMF = 40; both also: R024S485, contra
            case 22
                rabbit = 'R024';
                session = 'R024S402_TT1N2'; CF = 400; %y = [0 120]; BMF = 27; no contra
            case 23
                rabbit = 'R027';
                session = 'R027S019_TT3N1'; CF = 1000;
        end

        datafile = fullfile(base, fpath, rabbit, session);
        data_full = load(datafile, 'saved_data');
        data_full = data_full.saved_data;
        has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
            strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
            d.param.binmode == 2&&(d.param.spl==spl(spl_ind)||d.param.spl==spl2(spl_ind))&&mod(d.param.fpeak_mid, 200)==0 ...
            ,data_full);

        % Plot Spectral Centroid
        data = data_full(has_sc);
        num_DSIDs = length(data);
        for ind = 1:num_DSIDs
            DSID = data{ind}.param.dsid;

            [fpeaks,~,fpeaksi] = unique([data{ind}.param.list.fpeak].');
            num_fpeaks = length(fpeaks);
            dur = data{ind}.param.dur/1000; % stimulus duration in seconds.

            rate_size = [num_fpeaks,1]; % [num_F0s,num_Fps];
            spike_rates = data{ind}.spike_info.num_spikes_delayed/...
                (dur - data{ind}.spike_info.onsetWin/1000);

            if length(fpeaksi) == length(spike_rates)
                [rate,rate_std] = accumstats({fpeaksi},spike_rates, rate_size);
                rate = smooth(rate);

                if length(rate) == 40
                    rates(j,:) = [rate; 0];
                else
                    rates(j,:) = rate;
                end
                CFs(j) = CF;
                label(j) = {num2str(CF)};
                j = j+1;
            end
        end
    end

    [CFs,index] = sort(CFs);
    rates = rates(index,:);
    CFs_unique = unique(CFs);

     g = linspace(0, 0.9, length(CFs_unique))';
     r = g.^1.6;
     b = g.^2.1;
     mymap = [r g b];

    norm_fpeaks = linspace(200, 1800, 41);
    subplot(1, 3, spl_ind)
    hold on
    for iii = 1:length(CFs)
        x = find(CFs(iii)==CFs_unique);
        plot(norm_fpeaks,rates(iii,:),'LineWidth', 5, 'Color', mymap(x,:))
    end
    %plot(norm_fpeaks,rates([end-1, end],:),'LineWidth', 3, 'Color', 'b')
    legend(strcat(' ', num2str(CFs')), 'location', 'northeastoutside')
    if spl_ind == 2
        xlabel('Spectral Peak Frequency re. CF (Hz)')
    end
    xlim([200 1800]);
    xticks([200 400 600 800 1000 1200 1400 1600 1800])
    xticklabels({'-800', '', '-400', '', 'CF', '', '+400', '', '+800'})
    if spl_ind == 1
        ylabel('Avg. Rate (sp/s)')
    end
    ylim([0 110])
    %ylim([0 1])
    set(gca,'FontSize',20)
    %set(gca, 'XScale', 'log')
    box on
    grid on
    title([num2str(spl(spl_ind)) ' dB SPL'], 'FontSize', 22);
end

