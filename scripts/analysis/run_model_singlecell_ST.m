%% run_model_singlecell_ST.m
%
% This script runs and plots HSR, LSR, BE, and BS responses from the UR_EAR
% model for a single-cell response for synthetic timbre for different
% slopes and different F0s. 
%
% Other m-files required: generate_ST, modelAN, wrapperIC
%
% Author: J. Fritzinger
% Created: 2023-08-29; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------
clear


timerVal = tic;

figure('Position',[914,508,646,730]);
tiledlayout(5, 2, "TileSpacing","compact", "Padding","compact")

F0 = [100 200];
%spl = [43, 63, 83];
spl = 70;
g = [200 100 24];
for iF0 = 1:length(F0)
	for ispl = 1:length(g)

		%% Create stimulus

		% Stimulus parameters
		params.type = 'SPEC_slide';
		params.fpeak_mid = 1600;
		params.num_harms = 13;
		params.stp_otc = 41;
		params.Delta_F = F0(iF0);
		params.dur = 0.3;
		params.ramp_dur = 0.02;
		params.steps = 41;
		params.spl = spl;
		params.g = g(ispl);
		params.Fs = 100000;
		params.range = 3;
		params.version = 5;
		params.mnrep = 1;
		params.physio = 0;
		params.reptim = 0.6;

		% Generate sliding stimulus
		[params] = generate_ST(params);
		params.num_stim = size(params.stim, 1);


		%% Model
		CF = 1600;

		% Model parameters
		model_params.type = 'SFIE';
		model_params.range = 2; % 1 = population model, 2 = single cell model
		model_params.species = 1; % 1 = cat, 2 = human
		model_params.BMF = 63;
		model_params.CF_range = CF;
		model_params.num_CFs = 1;
		model_params.CFs = CF;
		model_params.nAN_fibers_per_CF = 5;
		model_params.cohc = 1; % (0-1 where 1 is normal)
		model_params.cihc = 1; % (0-1 where 1 is normal)
		model_params.nrep = 1; % how many times to run the AN model
		model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
		model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
		model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
		model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050

		% Run model
		model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
		AN_HSR = modelAN(params, model_params); % HSR for IC input
		model_params.fiberType = 1; % LSR
		AN_LSR = modelAN(params, model_params); % LSR for plotting
		SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output

		% % Run efferent model
		% nstim = size(params.stim,1);
		% CFs = CF;
		% nCFs = length(CFs);
		% SFIE.average_ic_sout_BE = zeros(nstim, nCFs);
		% AN_HSR.average_AN_sout = zeros(nstim, nCFs);
		% AN_LSR.average_AN_sout = zeros(nstim, nCFs);
		% for istim = 1:nstim
		% 	stimulus = params.stim(istim,:);
		% 
		% 	% Call model w/ efferent system enabled
		% 	[ihc_eff(istim,:,:), hsr_eff(istim,:,:), lsr_eff(istim,:,:), ic_eff(istim,:,:), gain_eff(istim,:,:)] = ...
			% 	sim_efferent_model(stimulus, CFs, moc_weight_wdr=2.0, moc_weight_ic=8.0, species=1);
		% 
		% 
		% 	% Average and store IC rate
		% 	%SFIE.average_ic_sout_BE(istim,:) = mean(ic_eff(istim,:,:),3);  % note: this includes onset!
		% 	AN_HSR.average_AN_sout(istim,:) = mean(hsr_eff(istim,:,:),3);
		% 	AN_LSR.average_AN_sout(istim,:) = mean(lsr_eff(istim,:,:),3);
		% end
		% SFIE = wrapperIC(hsr_eff, params, model_params); % SFIE output
		% %SFIE.ic_BE = ic_eff;

		%% Plot results

		xlimits = [params.fpeaks(1) params.fpeaks(end)];

		nexttile(iF0)
		[~, index] = min(abs([params.mlist.fpeak]-CF));
		y = fft(params.stim(index,:));
		m = abs(y);
		mdB = 20*log10(m);
		f = (0:length(y)-1)*params.Fs/length(y);
		hold on
		plot(f,mdB, 'LineWidth',2)
		xlim(xlimits);
		ylim([0, 80])
		%set(gca, 'XScale', 'log');
		grid on
		box on
		ylabel('Mag (dB SPL)')
		xlabel('Frequency (Hz)')
		title(['F0 = ' num2str(F0(iF0)) 'Hz'])


		nexttile(iF0+length(F0))
		hold on
		[avAN_HSR, stdAN_HSR] = plotST(params, AN_HSR.average_AN_sout, 0);
		errorbar(params.fpeaks, avAN_HSR, stdAN_HSR./sqrt(params.mnrep), 'LineWidth',2)
		%set(gca, 'XScale', 'log');
		ylabel('Rate (sp/sec)','fontsize',10)
		xlabel('CF (Hz)')
		title('AN HSR Response Profile','fontsize',10)
		xlim(xlimits);
		ylim([0 max(avAN_HSR)+5])
		box on
		grid on
% 		if ispl == length(spl)
% 			for iispl = 1:length(spl)
% 				leg{iispl} = [num2str(spl(iispl)) ' dB SPL'];
% 			end
% 			legend(leg, 'Location','best');
% 		end

		nexttile(iF0+length(F0)*2)
		hold on
		[avAN_LSR, stdAN_LSR] = plotST(params, AN_LSR.average_AN_sout, 0);
		%plotSyntheticTimbre_Windowed(params, SFIE.ic_BE, 1);
		errorbar(params.fpeaks, avAN_LSR, stdAN_LSR./sqrt(params.mnrep), 'LineWidth',2)
		ylim([0 max(avAN_LSR)+5])
		%set(gca, 'XScale', 'log');
		xlim(xlimits);
		ylabel('Rate (sp/sec)','fontsize',10)
		xlabel('Spectral Peak Frequency (Hz)')
		title('AN LSR Response Profile','fontsize',10)
		box on
		grid on

		nexttile(iF0+length(F0)*3)
		hold on
		[avIC_BE, stdIC_BE] = plotST(params, SFIE.average_ic_sout_BE, 0);
		%plotSyntheticTimbre_Windowed(params, SFIE.ic_BE, 1);
		errorbar(params.fpeaks, avIC_BE, stdIC_BE./sqrt(params.mnrep), 'LineWidth',2)
		ylim([0 max(avIC_BE)+5])
		xlim(xlimits);
		%set(gca, 'XScale', 'log');
		ylabel('Rate (sp/sec)','fontsize',10)
		xlabel('Spectral Peak Frequency (Hz)')
		title('BE Response Profile','fontsize',10)
		box on
		grid on

		nexttile(iF0+length(F0)*4)
		hold on
		[avIC_BS, stdIC_BS] = plotST(params, SFIE.average_ic_sout_BS, 0);
		%plotSyntheticTimbre_Windowed(params, SFIE.ic_BS, 1);
		%legend('0-100ms', '100-200ms', '200-300ms')
		
		errorbar(params.fpeaks, avIC_BS, stdIC_BS./sqrt(params.mnrep), 'LineWidth',2)
		xlim(xlimits);
		ylim([0 max(avIC_BS)+5])
		%set(gca, 'XScale', 'log');
		ylabel('Rate (sp/sec)','fontsize',10)
		xlabel('Spectral Peak Frequency (Hz)')
		title('BS Response Profile','fontsize',10)
		box on
		grid on

	end
	legend('200 db/oct', '100 db/oct', '24 db/oct');
end
elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])