function [DV, weighted_avg] = spec_centroid_model_mod(cond, stimulus, order)
% Plot AN and IC rates for N=# repetitions

Fs = 100e3; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):

% Simulation parameters
nAN_fibers_per_CF = 1; % Number of stimulus repetitions (only 1 rep is needed if probability of firing is displayed - more reps needed for looking at spike times)
num_CFs = 25;
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
CF_range = [125 3000]; % set a 'generic' default range of CFs
rng('default');     % added to accommodate older versions of Matlab, in cases for which rng has previously been run
rng('shuffle');     % seed the random number generator using time of day
CFs = logspace(log10(CF_range(1)),log10(CF_range(2)),num_CFs); % set range and resolution of CFs here

%% Run Model

for trial = 1:2  
    for inoise = 1:num_noises
        %% Set up and run the simulation
        dur = length(stimulus(trial,:))/Fs;
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
                case 2 %Monaural Simple Filter
                    ic_sout_BE = unitgain_bpFilter(an_sout,BMF,Fs);  % Now, call NEW unitgain BP filter to simulate bandpass IC cell with all BMF's
                    average_ic_sout_BE(inoise,iCF) = mean(ic_sout_BE(onset_num:end)); % averages the bandpass response over the stimulus duration
            end
        end % END OF CF LOOP
    end % noises
    
    %% Plot
%     switch cond % 1 = timbre, 2 = pitch
%         case 1
%             g = new_approach(average_ic_sout_BS, CFs);
%             plotModelResponses(trial, CFs, average_AN_sout(inoise,:), average_ic_sout_BE(inoise,:), ...
%                 average_ic_sout_BS(inoise,:), stimulus, Fs, Which_IC, cond, Fc, g)
%         case 2
%             plotModelResponses(trial, CFs, average_AN_sout(inoise,:), average_ic_sout_BE(inoise,:), ...
%                 average_ic_sout_BS(inoise,:), stimulus, Fs, Which_IC, cond, F0)
%            
%     end
    %% Calculate decision variable
    switch cond % 1 = timbre, 2 = pitch
        case 1    
            
            % DV 2
            %weighted_avg(trial) = (sum(CFs.*average_ic_sout_BS))/(sum(average_ic_sout_BS));
            
            % DV 3
%             y = supsmu(CFs, average_ic_sout_BS, 'Span', 0.5);
%             y_d = average_ic_sout_BS-y;
%             y_d(y_d<0)=0;
%             weighted_avg(order(trial)) = (sum(CFs.*y_d))/(sum(y_d));
            
            % new approach, Braden's DV 
            weighted_avg(trial) = new_approach(average_ic_sout_BS, CFs);
        case 2
            new = -average_ic_sout_BS(inoise,:);
            [~, indices] = findpeaks(new, 'MinPeakProminence', 2);
            BS_peaks = CFs(indices);
            below_Fc = find(BS_peaks<1000);
            if size(below_Fc', 1) == 1
                weighted_avg(trial) = BS_peaks(below_Fc(1));
            elseif size(below_Fc', 1) == 0
                weighted_avg(trial) = 0;
            else
                weighted_avg(trial) = BS_peaks(below_Fc(end))-BS_peaks(below_Fc(end-1));
            end
    end
end % conditions

%% Decision Variable
switch cond % 1 = timbre, 2 = pitch
    case 1
        DV = 1; % correct
        if order(1) == 1 % low played first 
            if weighted_avg(1) >= weighted_avg(2)
                DV = 0; % incorrect
            end
        elseif order(1) == 2 % high played first 
            if weighted_avg(2) >= weighted_avg(1)
                DV = 0; % incorrect
            end
        end
    case 2
        DV = 1;
        if weighted_avg(1) >= weighted_avg(2)
            DV = 0;
        end
end
end
