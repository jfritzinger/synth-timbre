% Model_Predictions_MTFs

%% Generate SFIE MTF
% J. Fritzinger
clear all 
close all

%% Analyze MTF data 
% Find MTF 
%rate_sm = smooth_rates(data{1,1}.rate,data{1,1}.rlb,data{1,1}.rub);
%[BMF,WMF,MTF_shape] = MTFclassification(rate_sm, data{1,1}.fms);

% Stimulus Parameters and Generation
%params = data{1,1}.param;
params.type = 'typMTFN';
params.ramp_dur = 0.05;
params.noise_state = 0;
params.noise_band = [100, 10000];
params.dur = 1; % s
params.reptim = 1.5;
params.fms = [2, 600, 3]; % fm_lo, fm_hi, steps per octave
params.mdepths = [0,0,1];
params.binmode = 2;
params.No = 30;
params.spl = 30;
params.all_mdepths = 0;
params.Fs = 100000;
params.nrep = 1;

% Model Parameters 
model_params.BMF = 100;
model_params.Fs = params.Fs;
model_params.species = 1;
model_params.CF = 1000;
% if strcmp(MTF_shape, 'BE')
%     model_params.BMF = BMF;
% elseif strcmp(MTF_shape, 'BS')
%     model_params.BMF = WMF;
% else
%     model_params.BMF = 100;
% end

% Generate stimulus 
if strcmp(params.type,  'typMTFN')
    [stim, ~, model_params.list, params] = generate_MTF_2(params);
    model_params.stim_list = 1:length([model_params.list.fm]);
    model_params.fpeaks = length([model_params.list.fm]);
    model_params.steps = size(model_params.list, 2);
end

% Run SFIE model 
[~, ~, avBE_SFIE, stdBE_SFIE, avBS_SFIE, stdBS_SFIE, stim_ref] = modelSingleCell_2(stim, model_params, params);

% Plot MTF
%data = data{1,1};
[~, ind] = sort([model_params.list.fm]);
figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,469,250,190]); % left bottom width height

%% Plot SFIE model 
avBE_SFIE = avBE_SFIE(ind); % sort
stdBE_SFIE = stdBE_SFIE(ind); % sort 
avBS_SFIE = avBS_SFIE(ind); % sort
stdBS_SFIE = stdBS_SFIE(ind); % sort 

hold on
line([1 params.all_fms(end)], [1 1]*avBE_SFIE(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
%errorbar(params.all_fms,avBE_SFIE, stdBE_SFIE,'.');
%line(params.all_fms,avBE_SFIE,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w');
plot(params.all_fms,smooth(avBE_SFIE),'-b', 'LineWidth', 2) % smoothed MTF
hold off
% Label the plots
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
xlabel('Modulation Freq (Hz)')
ylabel('Rate (sp/s)')
set(gca,'XTick',xtick,'XScale', 'log')
legend('Unmodulated', 'Location','northwest')
grid on
set(gca,'FontSize',20)
xticklabels([])
yticklabels([])
box on

figure;
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(gcf,'color','w');
set(gcf, 'position', [560,469,250,190]); % left bottom width height
hold on
line([1 params.all_fms(end)], [1 1]*avBS_SFIE(1),'Color',[0.7 0.7 0.7], 'LineWidth', 2);
%errorbar(params.all_fms,avBS_SFIE, stdBS_SFIE,'.');
%line(params.all_fms,avBS_SFIE,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w');
plot(params.all_fms,smooth(avBS_SFIE),'-b', 'LineWidth', 2) % smoothed MTF
hold off
% Label the plots
xtick = [1 2 5 10 20 50 100 200 500];
xlim(xtick([1 end]))
xlabel('Modulation Freq (Hz)')
ylabel('Rate (sp/s)')
ylim([0 34])
set(gca,'XTick',xtick,'XScale', 'log')
legend('Unmodulated', 'Location','northwest')
xticklabels([])
yticklabels([])
grid on
set(gca,'FontSize',20)
box on


%% Single-Cell Model 
function [avAN, stdAN, avBE, stdBE, avBS, stdBS, stim_ref] = modelSingleCell_2(stim, model_params, params)
% Returns the average rates and stds for AN, BE, and BS for a single CF
% J. Fritzinger, updated 2/1/2022

% Model Parameters
Fs = model_params.Fs; % Modeling sampling rate in Hz (must be 100, 200 or 500 kHz for AN model):
cohc =  1;          % (0-1 where 1 is normal)
cihc =  1;
fiberType =  3;     % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
implnt = 0;         % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
noiseType = 1;      % 0 for fixed fGn (1 for variable fGn) - this is the 'noise' associated with spontaneous activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
Which_IC = 1;       % 2 = ModFilt; 1 = SFIE model; BMF will be matched to talker's F0 below
numreps = 10;       %number of runs

% Run model
for i = 1
    BEICresp = zeros(numreps,model_params.steps);
    BSICresp = zeros(numreps,model_params.steps);
    ANresp = zeros(numreps,model_params.steps);
    for j = 1:numreps %for multiple runs
        
        if (strcmp(params.type, 'SPEC_slide') && strcmp(params.type, 'WB_noise'))...
            || strcmp(params.type, 'typMTFN') || strcmp(params.type,'Natural_Timbre')
            stim_ref = 0; % changed for WB noise code
        else
            [~, stim_index] = min(abs(Fc(i) - model_params.fpeaks));
            stim_ref(:,i) = stim(:,find(model_params.stim_list(i,:) == stim_index));
        end
        
        % Runs each presentation through AN model
        for n = 1:model_params.steps
            
            stimulus=(stim(:,n))';
            nrep = 1;%repitition of simulus in one run
            dur = length(stimulus)/Fs;             % duration of waveform in sec
            
            % Using ANModel_2014 (2-step process)
            vihc = model_IHC(stimulus, model_params.CF, nrep, 1/Fs, dur*1.2, cohc, cihc, model_params.species);
            num_AN_fibers = 10; %average across AN fibers -4/7/2020
            for iAN = 1:num_AN_fibers
                [an_sout,~,~] = model_Synapse(vihc,model_params.CF,nrep,1/Fs,fiberType,noiseType,implnt); % an_sout is the auditory-nerve synapse output - a rate vs. time function that could be used to drive a spike generator
                ANfibers_response(iAN,:) = an_sout;
            end
            
            an_sout = mean(ANfibers_response);%mean across independent AN fibers and reassign to an_sout - MA-4/1/20
            AN_average = mean(an_sout); % save mean rates for a plot of population AN response
            
            %IC model- using SFIE model
            [ic_sout_BE,ic_sout_BS,~] = SFIE_BE_BS_BMF(an_sout,model_params.BMF,Fs); %%%new model
            average_ic_sout_BE = mean(ic_sout_BE); % averages the bandpass response
            average_ic_sout_BS = mean(ic_sout_BS);
            
            % save IC & AN responses for each runs
            ANresp(j,model_params.stim_list(n))=AN_average;%AN response
            BEICresp(j,model_params.stim_list(n)) = average_ic_sout_BE;%IC band-enhanced cell response
            BSICresp(j,model_params.stim_list(n)) = average_ic_sout_BS;%IC band-suppressed cell response            
        end
    end
    
    %get average and standard deviation of results from multiple runs
    avAN(i,:)=mean(ANresp,1);
    stdAN(i,:)=std(ANresp,0,1);
    
    avBE(i,:)=mean(BEICresp,1);
    stdBE(i,:)=std(BEICresp,0,1);
    
    avBS(i,:)=mean(BSICresp,1);
    stdBS(i,:)=std(BSICresp,0,1);
    
end
end

%% Generate MTF stimulus
function [stimuli, stim_list, list, params] = generate_MTF_2(params)
% Modulation transfer function with band-limited noise carrier, for use with National Insturments hardware and naq program

% Calculate stimulus modulation depths
params.all_mdepths = params.mdepths(1):params.mdepths(3):params.mdepths(2);
nmdepths = length(params.all_mdepths);
fs = params.Fs;

% Calculate stimulus modulation frequencies
fm_lo = params.fms(1);
fm_hi = params.fms(2);
steps_per_octave = params.fms(3);
i = 1;
fm = fm_lo;
while fm < fm_hi
    all_fms(i) = fm;
    fm = fm * 2^(1/steps_per_octave); % next freq
    i = i+1;
end
if nmdepths == 1    % MTF at single mdepth
    all_fms_fig = [1.2 all_fms]; % Was a BUG: adding 2 Hz was trying to plot 0 Hz on log axis with the old setup
    params.all_fms = [0 all_fms];            % Add 0 Hz modulation frequency
end
nfms = length(params.all_fms);

% Estimate time required for this DSID
nstim = nfms * nmdepths;
time_required = nstim * params.nrep * params.reptim/60; %what is reptim
disp(['This will take ' num2str(time_required) ' min']);
npts = floor(params.dur*params.Fs);

% Create the Stimulus Gating function
npts = floor(params.dur*fs);
gate = tukeywin(npts,2*params.ramp_dur/params.dur); %raised cosine ramps
t = (0:1:npts-1)/fs;

% Construct filters
fn = fs/2;
b_bp = fir1(4000,params.noise_band/fn);         % 1. Bandpass filter
params.stim_filters.b_bp = b_bp;

% Begin calculating noise carrier
noise_spl = params.No+10*log10(fs/2); % SPL with desired spectrum level 
% use fs as bandwidth but will filtered later to get the correct overall level
noise_rms = 10^(noise_spl/20)*20e-6;                % RMS amplitude
if params.noise_state == 0                          % Frozen noise; random noise is handled in the stim loop
    rng(0);
    BBN_Pa = noise_rms*randn(npts,1);               % White noise (BW: 0 - fn) with appropriate spectrum level
    BBN_Pa_band = conv(BBN_Pa,b_bp,'same');         % Band-limited noise of appropriate spectrum level
    pre_mod_rms = rms(BBN_Pa_band);                 % RMS amplitude of noise waveform prior to modulation (target RMS)
end

% Generate stimuli for all presentations
stimuli = zeros(npts,nstim*params.nrep);
presentation = 0;

for irep = 1:params.nrep
    rng('shuffle');                               % For each rep, create random sequence of stimuli
    stim_list = randperm(nstim);
    for istim = 1:nstim
        presentation = presentation + 1;
        
        imdepth = mod(stim_list(istim),nmdepths)+1;             % "+1" starts values at 1 instead of 0
        mdepth = params.all_mdepths(imdepth);
        ifm = ceil((stim_list(istim))/nmdepths);
        fm = params.all_fms(ifm);
        
        % Finish calculating noise carrier waveform
        if params.noise_state == 0                              % Frozen noise carrier
            seed = 0;
        elseif params.noise_state == 1                          % Random noise
            seed=rng;
            BBN_Pa = noise_rms*randn(npts,1);                   % White noise (BW: 0 - fn) with appropriate spectrum level
            BBN_Pa_band = conv(BBN_Pa,b_bp,'same');             % Band-limited noise of appropriate spectrum level
            pre_mod_rms = rms(BBN_Pa_band);                     % RMS amplitude of noise waveform prior to modulation
        end
        
        % Calculate modulator tone
        if mdepth <= -33
            m = 0;                                              % Replace "-35 dB modulation" with an unmodulated stimulus
        else
            m = 10^(mdepth/20);                                 % Convert mdepth in dB into m
        end
        modulator = m*sin(2*pi*fm*t');
        
        % Modulate the noise carrier
        am_stim = (1 + modulator) .* BBN_Pa_band;
        post_mod_rms = rms(am_stim);
        stimulus = am_stim * pre_mod_rms / post_mod_rms;       % Scale RMS back to pre-mod RMS (i.e. the desired No)
		stimuli(:, presentation) = stimulus.*gate;  
        
        list(presentation).fm = fm;
        list(presentation).ifm = ifm;
        list(presentation).mdepth = mdepth;
        list(presentation).imdepth = imdepth;
        list(presentation).noise_seed = seed;
        list(presentation).rep = irep;
           
    end
end
end