%% Fig. 5: Robust to sound level
% Plots level differences, with a multi-unit physiology response
% J. Fritzinger, updated 2/23/2021
clear 

%% Run model

param = cell(4,1);
levels = [40 50 70 80];
for istim = 1:4
	param{istim}.fpeak_mid = 1200;
	param{istim}.Delta_F = 200;
	param{istim}.dur = 0.3;
	param{istim}.ramp_dur = 0.02;
	param{istim}.spl = levels(istim);
	param{istim}.g = 24;
	param{istim}.num_harms = 9;
	param{istim}.stp_otc = 41;
	param{istim}.Fs = 100000;
	param{istim}.mnrep = 1;
	param{istim}.physio = 0;
	param{istim} = generate_ST(param{istim});
	param{istim}.num_stim = size(param{istim}.stim, 1);
end

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

% Run model
SFIE = cell(4,1);
for istim = 1:4
	model_params.CF_range = param{istim}.fpeak_mid;
	model_params.CFs = param{istim}.fpeak_mid;

	AN = modelAN(param{istim}, model_params); % HSR for IC input
	SFIE{istim} = wrapperIC(AN.an_sout, param{istim}, model_params); % SFIE output

	if istim == 1
		[avAN_low, ~] = plotST(param{istim}, AN.average_AN_sout, 0);
		[avBS_low, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	elseif istim == 2
		[avAN_med, ~] = plotST(param{istim}, AN.average_AN_sout, 0);
		[avBS_med, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	elseif istim == 3
		[avAN_hi, ~] = plotST(param{istim}, AN.average_AN_sout, 0);
		[avBS_hi, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	else
		[avAN_xhi, ~] = plotST(param{istim}, AN.average_AN_sout, 0);
		[avBS_xhi, ~] = plotST(param{istim}, SFIE{istim}.average_ic_sout_BS, 0);
	end
end

%% Plots 

% Plot AN model
figure;
subplot(1, 3, 1)
hold on

plot(param{1}.fpeaks, avAN_low, 'LineWidth', 2, 'Color', '#99CEC5') 
plot(param{1}.fpeaks, avAN_med, 'LineWidth', 2, 'Color', '#44AA99')
plot(param{1}.fpeaks, avAN_hi, 'LineWidth', 2, 'Color', '#117733')
plot(param{1}.fpeaks, avAN_xhi, 'LineWidth', 2, 'Color', '#023A30') 

xticks([600 800 1000 1200 1400 1600 1800])
xlabel('Frequency (Hz)')
xlim([600 1800])
grid on
box on
title('AN Model Response')
ylabel('Rate (sp/sec)')
ylim([0 250])
legend('40 dB SPL','50 dB SPL','70 dB SPL', '80 dB SPL', 'location', 'southwest')
set(gca,'fontsize',14)

% Plot BS model
subplot(1, 3, 2)
hold on

plot(param{1}.fpeaks, avBS_low, 'LineWidth', 2, 'Color', '#95D3F7') 
plot(param{1}.fpeaks, avBS_med, 'LineWidth', 2, 'Color', '#2A8ABF')
plot(param{1}.fpeaks, avBS_hi, 'LineWidth', 2, 'Color', '#2F66A9')
plot(param{1}.fpeaks, avBS_xhi, 'LineWidth', 2, 'Color', '#140561') 

xticks([600 800 1000 1200 1400 1600 1800])
xlabel('Frequency (Hz)')
xlim([600 1800])
grid on
box on
title('IC BS Model Response')
ylabel('Rate (sp/sec)')
ylim([0 50])
legend('40 dB SPL','50 dB SPL','70 dB SPL', '80 dB SPL', 'location', 'southwest')
set(gca,'fontsize',14)

% IC Analysis
clear

base = getPaths();
fpath = 'data/aro-2021/Figure 6/R024S313_TT3N1_ds';

% Load
for ii = 1:4
    data(ii) = load(fullfile(base, [fpath num2str(7+ii) '.mat'])); 
end

for ii = 1:4
    label(ii, :) = [num2str(data(ii).param.spl) ' dB SPL'];
    
    [fpeaks,~,fpeaksi] = unique([data(ii).param.list.fpeak].');
    num_fpeaks = length(fpeaks);
    dur = data(ii).param.dur/1000; % stimulus duration in seconds.
    
    rate_size = [num_fpeaks,1];
    spike_rates = data(ii).spike_info.num_spikes_delayed/...
        (dur - data(ii).spike_info.onsetWin/1000);
    [rate(ii,:),~] = accumstats({fpeaksi},spike_rates, rate_size);    
end

% Plot IC responses
subplot(1, 3, 3)
hold on
plot(fpeaks,rate(2,:), 'LineWidth', 2, 'DisplayName', label(2,:), 'color', '#ECC262');
plot(fpeaks,rate(3,:), 'LineWidth', 2, 'DisplayName', label(3,:), 'color', '#E68940');
plot(fpeaks,rate(1,:), 'LineWidth', 2, 'DisplayName', label(1,:), 'color', '#BD5503');
plot(fpeaks,rate(4,:), 'LineWidth', 2, 'DisplayName', label(4,:), 'color', '#7D3700');

xticks([600 800 1000 1200 1400 1600 1800])
xlabel('Frequency (Hz)')
xlim([600 1800])
grid on
box on
title('Single-Unit IC Response')
ylabel('Rate (sp/sec)')
ylim([0 200])
set(gca,'fontsize',14)
hold off
set(gcf,'color','w');
legend('location', 'southwest')
