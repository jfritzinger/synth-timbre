%% Figure 3: 3D Flucuations
% J. Fritzinger, updated 2/23/2021
clear

%% Parameters

% Model parameters
Fs = 100000;
nrep = 30; % repetitions averaged for profile plots
species = 2; % human
implnt = 0;         % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
noiseType = 1;      % 0 for fixed fGn (1 for variable fGn) - this is the 'noise' associated with spontaneous activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
fiberType = 3;      % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
cohc = 1;
cihc = 1;

% Stimulus parameters
param.fpeak_mid = 1200;
param.Delta_F = 200;
param.dur = 0.3;
param.ramp_dur = 0.02;
param.spl = 70;
param.g = 24;
param.num_harms = 11;
param.stp_otc = 1;
param.Fs = 100000;
param.mnrep = 1;
param.physio = 0;
param = generate_ST(param);
stimulus = param.stim;
BMF = 100;
onset_num = 1;

% Figure parameters
CF_list = [800 1000 1200 1400 1600];
start_time = 0.1; % s - for plotting
start  = start_time * Fs;
stop_time = 0.15; % s
stop = stop_time * Fs;
b_LP = fir1(5000,100/(Fs/2),'low'); % LP at 100 Hz, for envelope plots


%% Model and Plotting
iplot = 0; % step thru CFs within each panel
inoise = 1;
gamma_param.srate = Fs;
RMS_DV = zeros(nrep, 1);
env_slope_DV = zeros(nrep, 1);
rate_DV = zeros(nrep, 1);
avBE = zeros(1, length(CF_list));
avBS = zeros(1, length(CF_list));
fluctuation_amp = zeros(length(CF_list),1);
mean_rate = zeros(length(CF_list),1);
mean_RMS = zeros(length(CF_list),1);
for CF_plot = CF_list  % for TIN response to 1000 Hz tone
    iplot = iplot + 1;
    for irep = 1:nrep % do nrep repetitions - plot the last one, but use all 10 for averages
        impaired = 0; % not impaired!

		gamma_param.fc = CF_plot; % CF
		excit = gamma_filt(stimulus,gamma_param,impaired);
        RMS_DV(irep) = rms(excit); %mean(an_sout(start:stop))
        
        % AN model responses
        vihc = model_IHC(stimulus,CF_plot,nrep,1/Fs,param.dur*1.2,cohc,cihc,species);
        [an_sout,~,~] = model_Synapse(vihc,CF_plot,nrep,1/Fs,fiberType,noiseType,implnt); % an_sout is the auditory-nerve synapse output - a rate vs. time function that could be used to drive a spike generator
        an_sout_LP = conv(an_sout,b_LP,'Same'); % smooth AN PSTH to remove much of fine-structure, so it doesn't dominate a simple "enevelope slope" DV
        env_an_sout = envelope(an_sout_LP);
        env_slope_DV(irep) = mean(abs(diff(env_an_sout)*Fs))/mean(an_sout); %mean(abs(diff(env_an_sout(start:stop))*Fs))/mean(an_sout(start:stop));
        rate_DV(irep) = mean(an_sout); % mean(an_sout(start:stop));
        % RMS_DV(irep) = mean(an_sout(start:stop));
        
    end % last rep will be plotted below, but profiles will be based on means
    
    %% Make 3D plot
    
    % IC model
    [ic_sout_BE,ic_sout_BS,cn_sout_contra] = SFIE_BE_BS_BMF(an_sout, BMF, Fs);
    avBE(inoise,iplot) = mean(ic_sout_BE(onset_num:end));      % averages the bandpass response
    avBS(inoise,iplot) = mean(ic_sout_BS(onset_num:end));
    
    t = (1:length(an_sout))/Fs; % time vector for plots
    
    figure(12)  % AN responses  
    plot3(t(start:stop)*1e3,(CF_plot)*ones(1,length(t(start:stop))),an_sout(start:stop),'Color', '#44AA99','linewidth',1.5)
    hold on
    set(gca,'fontsize',14,'ytick',CF_list,'xtick',[100 125 150],'Ydir','reverse')
    
    % Duplicate the AN responses for Fluctuation profile plot
    plot3(t(start:stop)*1e3,(CF_plot)*ones(1,length(t(start:stop))),an_sout(start:stop),'Color', '#44AA99','linewidth',1)
    [peaks_vals,peak_locations] = findpeaks(an_sout(start:stop),Fs,'MinPeakProminence',100);
    plot3((peak_locations + start_time)*1e3,(CF_plot)*ones(1,length(peak_locations)),peaks_vals,'Color', '#117733','linewidth',4) % add envelope
    
    xlim([0.8*start_time*1e3, 1.2*stop_time*1e3])
    xlim([start_time*1e3, 1.2*stop_time*1e3])
    zlim([0 1500]); % rate
    %ylim([(0.9 * min(CF_list)) (1.1*max(CF_list))]); % freq
    ylim([(min(CF_list)) (max(CF_list))]); % freq
    
    % compute profiles by averaging across reps
    fluctuation_amp(iplot) = mean(env_slope_DV); %max(peaks_vals) - min(peaks_vals);
    mean_rate(iplot) = mean(rate_DV);
    mean_RMS(iplot) = mean(RMS_DV);
    
end

%% Add fluctuation profile plot
figure(12)
fluct_plot_position = 1.1*stop_time*1e3;
fluct_scalar = 25;
rate_scalar = 8;

% Plots the IC average rate profile
%plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list, fluct_scalar * fluctuation_amp,'o-','linewidth',3)
plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list, fluct_scalar * avBE,'o-','linewidth',3, 'Color', [0, 0.4470, 0.7410])
% plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list, rate_scalar * (mean_rate - 100),'go-','linewidth',3)

% Plots the dotted line bars
plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list,zeros(1,length(CF_list)),'k:','linewidth',2)
%plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list,1200*ones(1,length(CF_list)),'k:','linewidth',2); % plot a standard max line to help compare across panels
%plot3(fluct_plot_position*ones(1,length(CF_list)), CF_list,max(fluctuation_amp)*ones(1,length(CF_list)),'k:','linewidth',2)

% Plots the sticks
for ii = 1:length(CF_list)
    plot3([fluct_plot_position, fluct_plot_position],[CF_list(ii) CF_list(ii)],[0,fluct_scalar * avBE(ii)],'k','linewidth',1)
end
set(gca,'linewidth',2)
view([-50 65]) %([-55, 60])

%% Plot the stimulus
% stim_plot_position = 0.8*start_time*1e3;
% stim_scalar = 15;
% CF_list_harms = [800 1000 1200 1400 1600];
% 
% % Stimulus Creation
% harmonics = F0:F0:10000; % component freqs for the central stimulus, when this_fpeak = CF
% num_harmonics = length(harmonics);
% npts = dur * Fs; % # pts in stimulus
% t = (0:(npts-1))/Fs; % time vector
% component_scales_linear = 10.^(-1*abs(log2(harmonics/Fc)*G)/20); % one set of scales for the center triangle, i.e. when this_fpeak = CF
% stimulus = zeros(1,npts);
% for iharm = 1:num_harmonics
%     comp_freq = harmonics(iharm);
%     component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
%     stimulus = stimulus + component;          %Add component to interval
% end
% Level_scale = 20e-6*10.^(stimdB/20) * (1/rms(stimulus)); % overall linear scalar to bring this centered stimulus up to stimdB
% component_scales_linear = Level_scale * component_scales_linear; % include dB scaling into the set of harmonic component scalars
% 
% % Now, make the stimulus for this_fpeak
% ii = 1;
% for iharm = 1:num_harmonics
%     stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
%     
%     % Plot each stimulus
%     y = fft(stim);
%     mdB = 20*log10(abs(y));
%     level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);    
%     if harmonics(iharm) == CF_list_harms(ii) && ii < length(CF_list_harms)
%         plot3([stim_plot_position, stim_plot_position],[CF_list_harms(ii) CF_list_harms(ii)],...
%         [0,stim_scalar * level(ii)], 'Color', '#882255','linewidth',4)
%         ii = ii + 1;
%     end  
% end
% Add frequency lines 
for ii = 1:length(CF_list)
    plot3([start_time, fluct_plot_position],[CF_list(ii) CF_list(ii)],[0, 0],'k','linewidth',1)
end 

% Add title
set(gca,'fontsize',14)
set(gcf,'color','w');
%title('Spectral Centroid')
%xlabel('Time (ms)')
xticks([100 150])
xticklabels([0 50])
%zlabel('Rate (sp/s)')
zlim([0 1200])
%ylabel('Frequency (Hz)')

