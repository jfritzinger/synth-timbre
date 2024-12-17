%% run_latmodel_RVF.m
clear

%% Generate stimulus 

% RVF test 
params.type = 'RVF';
params.spl = 60; %dB SPL - changed to 60 dB SPL on 3/29/23 (PM)
params.binmode = 2; %0 - ipsi, 1 - contra, 2 - binaural
params.vels = [0.5:0.25:2,2.5:0.5:3,4:9]; % Velocities (kHz/ms)
params.per_reps = 1; % Period Reps (# consecutive chirp repeats)
params.seg_reps = 8; % Segment Reps (# presentations of each F0/C within a stimulus rep
params.nrep = 1; % Stimulus Reps (# identical stimulus presentations)
params.Fs = 100000;
params = generate_RVF(params);
params.num_stim = size(params.stim, 1);

% % SCHR
% params.type = 'SCHR'; % was "LNN" and before that "shiftTIN_CB"
% params.subtype = 'Regular'; % 'Sparse3' 'UpDown' << set below
% params.version = 1;
% params.dur = 0.4; % in the paper, some experiments use longer maskers than tones
% params.reptim = 1.0;
% params.ramp_dur = 0.025;
% params.binmode = 2; % binaural
% params.stimdB = 60;
% params.F0 = [50 100 200 400];
% params.C = -1.0:0.5:1.0;
% params.mnrep = 5;
% params.nrep = 5;
% params.Fs = 100000;
% [params] = generate_SCHR(params);
% params.num_stim = size(params.stim, 1);


%% Run model 

% Parameters for the lateral inhibition strength, delay, and CF range 
paramS1 = 0.5; % Strength of off-CF inhibition (recommend 0.3-0.5)
paramS2 = 0.0; % Strength of off-CF inhibition (recommend 0.3-0.5)
paramCF = 0.75; % CF of off-CF inhibitory pathways, octaves (recommend 0.5-1)
paramD = 0.002; % Off-CF inhibition delay, ms (has little effect on responses)
CF = 1500;

% AN model parameters 
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.nAN_fibers_per_CF = 1;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 1; % how many times to run the AN model
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.CFs = [CF*2^(-1*paramCF), CF, CF*2^paramCF];
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050

% IC model parameters 
model_params.type = 'Lateral Model';
model_params.config_type = 'BS inhibited by off-CF BS';
model_params.BMF = 100;
lm_params = [paramS1 paramS2 paramD];

% Runs AN lateral code 
AN = modelLateralAN(params, model_params);

% Runs lateral model SFIE code 
lateral_model = modelLateralSFIE(params, model_params, AN.an_sout, ...
	AN.an_sout_lo, AN.an_sout_hi,'CS_params', lm_params);

%% Plot model response 

% RVF
plotRVF(params, lateral_model.ic)

% % SCHR
% plotSCHR(params, lateral_model.ic)


