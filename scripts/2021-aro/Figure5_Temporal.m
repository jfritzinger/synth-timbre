%% Fig. 5 AN and IC temporal and average rate plots 
% This was going to be a temporal analysis figure, but was not finished. 
% J. Fritzinger, updated 2/23/2021
clear 

%% Model

% Stimulus parameters
param.fpeak_mid = 1100;
param.Delta_F = 200;
param.dur = 0.3;
param.ramp_dur = 0.02;
param.spl = 70;
param.g = 24;
param.num_harms = 9;
param.stp_otc = 41;
param.Fs = 100000;
param.mnrep = 1;
param.physio = 0;
param = generate_ST(param);
param.num_stim = size(param.stim, 1);

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 3; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.CF_range = param.fpeak_mid;
model_params.CFs = param.fpeak_mid;

% Run model
AN = modelAN(param, model_params); % HSR for IC input
SFIE = wrapperIC(AN.an_sout, param, model_params); % SFIE output
SFIE.ic_BS = squeeze(SFIE.ic_BS);

[avAN, ~] = plotST(param, AN.average_AN_sout, 0);
[avBS, ~] = plotST(param, SFIE.average_ic_sout_BS, 0);

%% Sync coefficient / vector strength
% Using Su & Delgutte 2019 

% Plot single fiber - time domain
figure;
t = linspace(0, length(SFIE.ic_BS)/param.Fs, length(SFIE.ic_BS));
plot(t, SFIE.ic_BS(10,:));
hold on
plot(t, SFIE.ic_BS(20,:));
plot(t, SFIE.ic_BS(15,:));
plot(t, SFIE.ic_BS(8,:));
ylabel('Spikes?') % don't know 
xlabel('Times (s)') % shouldn't this be 0.3s? 


% Period histogram (F0, envelope of stimulus)

% theta_i =  % phase of the ith spike relative to F0


% Compute vector strength to F0
% sum_cos = sum(cos(theta));  % 1xsteps
% sum_sin = sum(sin(theta));  
% VS = 1/avBS * sqrt(sum_cos.^2+sum_sin.^2); % vector strength 


% Plot VS to F0 vs. CF/Fc
% figure;
% plot(fpeaks, VS)

%% Plots 
figure(5);
set(gcf,'color','w');
xlimits = [min(param.fpeaks) max(param.fpeaks)];

% Top Left: AN, spectrogram
subplot(2, 3, 1)
title('Spectrogram')
grid on
box on

%[sg,Ftmp,Ttmp] = spectrogram(icBS(1,:),hanning(1000),[],2000,Fs); % 50% overlap (using even #'s to avoid "beating" with F0 for speech)
%spec_image = pcolor(Ttmp,Ftmp,abs(sg));

% Top Middle: AN, temporal
subplot(2, 3, 2)
title('Temporal')
grid on
box on

% Top Right: AN, average rate 
subplot(2, 3, 3)
plot(param.fpeaks, avAN, 'LineWidth', 2, 'Color', '#44AA99');
title('Average Rate')
grid on
box on
xlim(xlimits)
ylim([0 250])

% Bottom Left: IC, temporal
subplot(2, 3, 5)
xlabel('Frequency (Hz)')
grid on
box on

% Bottom Right: IC, average rate 
subplot(2, 3, 6)
plot(param.fpeaks, avBS, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
xlabel('Frequency (Hz)')
grid on
box on
xlim(xlimits)
ylim([0 50])

