function [DV, calc_DVs] = runModel(cond, stimulus, Fc, F0, Fc_order, F0_order, plot)
% Plot AN and IC rates for N=# repetitions
% J. Fritzinger, updated 1/6/2021
%
% Inputs: cond: pitch or timbre
%         stimulus: two stimuli 
%         Fc: spectral centroids
%         F0: fundamental frequencies
%         Fc_order: spectral centroid order
%         F0_order: fundamental freq order
%         plot: 0 for no plot, 1 for plotting


Fs = 100e3; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):

% Simulation parameters
nAN_fibers_per_CF = 5; % Number of stimulus repetitions (only 1 rep is needed if probability of firing is displayed - more reps needed for looking at spike times)
num_CFs = 50;
num_noises = 1;         % # of different noises, results will be averaged across them

%% Model Parameters

% Model parameters
cohc =  1;          % (0-1 where 1 is normal)
cihc =  1;
nrep = 1;
fiberType =  3;     % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
implnt = 0;         % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
noiseType = 1;      % 0 for fixed fGn (1 for variable fGn) - this is the 'noise' associated with spontaneous activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
species =  2;       % 1=cat; 2=human AN model parameters (with Shera tuning sharpness)
Which_IC = 1;       % 2 = ModFilt; 1 = SFIE model; BMF will be matched to talker's F0 below
BMF = 100;

% CF parameters
CF_range = [125 5000]; % set a 'generic' default range of CFs
rng('default');     % added to accommodate older versions of Matlab, in cases for which rng has previously been run
rng('shuffle');     % seed the random number generator using time of day
CFs = logspace(log10(CF_range(1)),log10(CF_range(2)),num_CFs); % set range and resolution of CFs here

%% Run Model

for trial = 1:2   
    for inoise = 1:num_noises        
        %% Set up and run the simulation
        dur = length(stimulus(trial, :))/Fs;
        onset_num = 1; %round(0.010*Fs); %1;  % 1st point that will be included in analyzed response  (allows exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050*Fs;)
        
        % Loop through CFs
        for iCF = 1:length(CFs)
            CF = CFs(iCF); % CF in Hz;
            vihc = model_IHC(stimulus(trial,:),CF,nrep,1/Fs,dur*1.2,cohc,cihc,species);
            for ifiber = 1:nAN_fibers_per_CF
                [an_sout_one_fiber(ifiber,:),~,~] = model_Synapse(vihc,CF,nrep,1/Fs,fiberType,noiseType,implnt); % an_sout is the auditory-nerve synapse output - a rate vs. time function that could be used to drive a spike generator
            end
            if nAN_fibers_per_CF > 1
                an_sout = mean(an_sout_one_fiber); % mean across n independent AN fibers with same CF (!! Note: if thre's only one AN fiber, this takes the wrong mean!)
            else
                an_sout = an_sout_one_fiber;
            end
            average_AN_sout(inoise,iCF) = mean(an_sout(onset_num:end)); % save mean rates for a plot of population AN response
            
            switch Which_IC
                case 1 %Monaural SFIE
                    [ic_sout_BE,ic_sout_BS,cn_sout_contra] = SFIE_BE_BS_BMF(an_sout, BMF, Fs);
                    average_ic_sout_BE(inoise,iCF) = mean(ic_sout_BE(onset_num:end));      % averages the bandpass response
                    average_ic_sout_BS(inoise,iCF) = mean(ic_sout_BS(onset_num:end));
                    all_ic_sout_BS(iCF, :) = ic_sout_BS(onset_num:end);
                case 2 %Monaural Simple Filter
                    ic_sout_BE = unitgain_bpFilter(an_sout,BMF,Fs);  % Now, call NEW unitgain BP filter to simulate bandpass IC cell with all BMF's
                    average_ic_sout_BE(inoise,iCF) = mean(ic_sout_BE(onset_num:end)); % averages the bandpass response over the stimulus duration
            end
        end % END OF CF LOOP
    end % noises
    
    %% Plot
    if plot == 1
        plotModelResponses(trial, CFs, average_AN_sout(inoise,:), average_ic_sout_BE(inoise,:), ...
            average_ic_sout_BS(inoise,:), stimulus, Fs, Which_IC, cond, Fc, F0, Fc_order, F0_order)
    end 
    
    %% Calculate decision variable
    DV_type = 5;
    calc_DV = generateDV(cond, DV_type, CFs, average_ic_sout_BS, all_ic_sout_BS, Fs);
    % fix bug so that calc_DVs are in order  
    switch cond
        case 1
            if Fc_order(1)  == 1 % low played first
                calc_DVs(trial) = calc_DV;
            elseif Fc_order(1)  == 2 % high played first
                if trial == 1
                    calc_DVs(2) = calc_DV;
                else
                    calc_DVs(1) = calc_DV;
                end
            end
        case 2
            if F0_order(1)  == 1 % low played first
                calc_DVs(trial) = calc_DV;
            elseif F0_order(1)  == 2 % high played first
                if trial == 1
                    calc_DVs(2) = calc_DV;
                else
                    calc_DVs(1) = calc_DV;
                end
            end
    end
end % conditions

%% Decision Variable
switch cond % 1 = timbre, 2 = pitch
    case 1
        DV = 1; % correct
            if calc_DVs(1) >= calc_DVs(2)
                DV = 0; % incorrect
            end
    case 2
        DV = 1;
        if calc_DVs(1) >= calc_DVs(2)
            DV = 0;
        end
end
end
