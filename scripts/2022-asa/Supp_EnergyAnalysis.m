% Supp_EnergyAnalysis 
clear all 
clear all 

%% Load data
rabbit = 'R024';
session = 'R024S470_TT2N1.mat';
if ismac
    fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Physiology/01 - Data/';
    slash = '/';
else
    fpath = 'C:\Users\johan\Box\02 - Physiology\01 - Data\';
    fpath = 'C:\Users\jfritzinger\Box\02 - Physiology\01 - Data\';
    slash = '\';
end

datafile = [fpath rabbit slash session];
data_full = load(datafile, 'saved_data');
data_full = data_full.saved_data;
has_sc = cellfun(@(d)strcmp(d.param.type,'SPEC_slide')&&...
    strcmp(d.param.SPEC_slide_type,'Spectral_Centroid')&&...
    d.param.spl ==40,data_full);
%load Fig3_BSModel
%load Fig3_Data
CF = 1150;
CFs = 1150;
BMF = 100;

%% 

data = data_full(has_sc);
num_DSIDs = length(data);
label_ind = 1;
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

        % Plot Data
        hold on
        h = errorbar(fpeaks,rate,rate_std/(sqrt(data{ind}.param.nrep)),'LineWidth', 1.5);
        label(label_ind) = {[num2str(data{ind}.param.spl) ' dB SPL']};
        label_ind = label_ind+1;
    end
end

%% Create stimulus

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
[params.freq_lo, params.freq_hi] = get_freq_limits(params.Fc, params.num_harms, params.F0);
if data{ind}.param.spl < 100
    params.stimdB = data{ind}.param.spl;
else
    params.stimdB = 100;
end

[stim, model_params.stim_list, model_params.fpeaks] = generate_spectralcentroid(params);


%% Run SFIE Model 
% 
% model_params.species = 1;
% model_params.steps = params.steps;
% model_params.CF = CF;
% model_params.BMF = BMF;
% [~, ~, avBE, stdBE, avBS, stdBS, ~] = modelSingleCell(stim, model_params, params);
% 
% for ind = 1:num_DSIDs
%     R_int = corrcoef(rate(:, ind),avBS(ind,:)).^2;
%     R(ind) = R_int(1, 2);
% 
%     disp(['Variance explained = ' num2str(R(ind))]);
% end

%% Energy Filter
Fs = params.Fs;
gamma_param.srate = Fs;
tvals = (1:length(stim))/Fs;
gamma_IF_reg = zeros(length(CFs),length(tvals));
impaired = 0; % 0 = not impaired; 1 = 'impaired'
for istim = 1:size(stim, 2) 
    gamma_param.fc = CFs; % CF
    pin_gamma(:,istim) = gamma_filt(stim(:,istim),gamma_param,impaired);
    [~,peak_indices_gamma] = findpeaks(diff(pin_gamma(:,istim))); % Find intervals between zero-crossings
    %    gamma_IF_irreg = (1./diff(peak_indices_gamma/Fs)) - CFs(iCF); % save IF RELATIVE to CF << for plotting
    gamma_IF_irreg = 1./diff(peak_indices_gamma/Fs); % save IF RELATIVE to CF
    %%%%%%%%%%subplot(4,1,4); plot(peak_indices_gamma(2:end)/Fs,gamma_IF,'r'); xlim([0 dur])
    irreg_vals = peak_indices_gamma(2:end)/Fs; % irregularly spaced time values
    gamma_IF_reg = interp1(irreg_vals,gamma_IF_irreg,tvals,'spline'); % Interpolate to get equally spaced points
    gamma_Fm_reg = Fs * diff( gamma_IF_reg); % Convert IF estimate to Fm estimate (dF/dt)
    gamma_Fm_reg_norm = Fs * diff( gamma_IF_reg)/CFs; % dF/dt, normalized by CF
    
    % Analysis
    y2 = fft(stim(:,istim));
    m = abs(y2);
    mdB = 20*log10(m);
    f = (0:length(y2)-1)*Fs/length(y2);
    mdB(mdB<0) = 0;
    f(f>Fs/2) = [];
    mdB = mdB(1:length(f));
    
%     % Plot 
%     figure;
%     subplot(1, 2, 1)
%     semilogx(f,mdB);
%     xlim([50 5000])
%     ylim([0 70])
%     grid on
%     set(gca, 'XScale', 'log')
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude (dB SPL)')
%     xticks([100 200 500 1000 2000 5000 10000 20000])
%     box on
    
    % Analysis
    y2 = fft(pin_gamma(:,istim));
    m = abs(y2);
    mdB = 20*log10(m);
    f = (0:length(y2)-1)*Fs/length(y2);
    mdB(mdB<0) = 0;
    f(f>Fs/2) = [];
    mag_dB(:,istim) = mdB(1:length(f));
    
%      % Plot
%     subplot(1, 2, 2)
%     semilogx(f,mag_dB(:, istim));
%     xlim([200 7000])
%     ylim([0 70])
%     grid on
%     set(gca, 'XScale', 'log')
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude (dB SPL)')
%     xticks([200 500 1000 2000 5000 10000 20000])
%     box on
%     title(['Gamma Filtered ' ])
    
end


%% Energy Calculation

for ii = 1:size(stim, 2)
    energy_filt(ii) = sum(abs(pin_gamma(:,ii)).^2);
end

plot(fpeaks, energy_filt);

% for ii = 1:size(stimulus, 1)
%     %energy_filt(ii, 1) = sum(abs(stimulus(ii,:)).^2);
%     energy_filt(ii) = sum(abs(mag_dB(ii,:)).^2);
% end

R = corrcoef(rate,energy_filt).^2;
R = R(1, 2);

disp(['Variance explained by energy = ' num2str(R)]);