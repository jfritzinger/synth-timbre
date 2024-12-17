%% NaturalTimbre_ModelPredictions 
% J. Fritzinger, updated 9/26/23
%
% This script loads in an example natural timbre response and plots model
% predictions (including variance explained) 


%% Load in info


if ismac
    stim_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/03 - Physiology/07 - Natural Timbre/Cut Waveforms/';
else
    stim_dir = 'C:\Users\jfritzinger\Box\03 - Physiology\07 - Natural Timbre\Cut Waveforms';
end

% Load in tuning sheet
if ismac
    tuning_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/03 - Physiology/07 - Natural Timbre/Cut Waveforms/';
else
    tuning_dir = 'C:\Users\jfritzinger\Box\03 - Physiology\07 - Natural Timbre\Cut Waveforms\';
end
tuning = readtable([tuning_dir 'Tuning.xlsx']);

%% Load in Data

% Identify current machine, which changes paths
[userid, base_dir, fpath, report_path] = findPaths();

% Reads in stimuli_sessions
if ismac
	filename = '/Volumes/CarneyLab/Rabbit_data/session_table.xlsx';
else
	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
end
celldata = readtable(filename, 'PreserveVariableNames',true);

% BE Example, R027S047
stimuli_of_interest = {'Natural_Timbre'};
rabbit = 25;
sesh = 595;
TT = 3;
N = 1;
CF = 1760;
MTF_shape = 'BS';

interest = celldata.Rabbit == rabbit & celldata.Session == sesh & ...
	celldata.Tetrode == TT & celldata.Neuron == N;
interest_ind = find(interest);
num_clusters = sum(interest);
session = sprintf('R%03dS%03d', celldata.Rabbit(interest_ind), celldata.Session(interest_ind));
rab_num = num2str(celldata.Rabbit(interest_ind));
if ispc
	rab_str =  ['R0' num2str(rabbit)];
	base_dir = base_dir{contains(base_dir,rab_str)};
	session_dir = fullfile(base_dir, session);
else
	session_dir = fullfile(base_dir, ['R0' num2str(rabbit)], session);
end

% Load data
[clusters, params, stims] = loadPhysiologySession(session_dir, session, userid);
cluster = clusters([clusters.tetrode] == TT); % Select tetrode
cluster = cluster([cluster.neuron] == N); % Select neuron

% Find binaural/contra datasets
is_contra = cellfun(@(p)isfield(p,'binmode')&&p.binmode ...
	== 1,params); % Sort data into binaural and contra
is_ipsi = cellfun(@(p)isfield(p,'binmode')&&p.binmode ...
	== 0,params); % Sort data ipsi
binmodes = {'Binaural', 'Contra', 'Ipsi'};
binmode = 'Binaural';
BIN = find(ismember(binmodes, binmode));
params_bin = params(~is_contra);

% Get which datasets have which sessions
has_nattimbre = cellfun(@(p)strcmp(p.type,'Natural_Timbre'),params_bin);

%% Physiology Analysis

num_DSIDs = sum(has_nattimbre);
data_NT = cell(num_DSIDs, 1);
params_NT = params_bin(has_nattimbre);
DSIDs = find(has_nattimbre);
for ind = 1 %1:num_DSIDst

	param = params_bin(DSIDs(ind));
	data = [];
	
	% Analysis
	[param, fig, data] = plotPhysNatTimbre(cluster, param, stims, 0, data);

end



%% Model - SFIE 

% Stimulus parameters
% params.type = 'Natural_Timbre';
% params.Fs = 100000;
% params.version = 2; % JBF 7/12/2021, added param for target input
% params.signal_spls = 70;
% params.noise_spls = 0;
% params.nrep = 5;
% params.noise_shape = [];
% params.noise_state = [];
% params.signal_ramp_dur = 0.02;
% params.signal_onset_delay = 0;
% params.binmode = 2;
% params.noise_band = [100 8000];             % Changed from [50 10000] on 12/7/15
% params.noise_ramp_dur = 0.02;
% params.target = 'Bassoon';
% params.reptim = 0.6;

% Parameters
params_NT{2}.Fs = 100000;
params_NT{2}.stim_dir = stim_dir;
params_NT{2}.mnrep = 1;

% Generate stimuli
[params_NT{2}] = generate_naturaltimbre(params_NT{2});
params_NT{2}.num_stim = size(params_NT{2}.stim, 1);

% Model parameters
model_params.BMF = 100;
model_params.species = 1;
model_params.CF = CF;
model_params.CF_range = CF;
model_params.num_CFs = 1;
model_params.CFs = CF;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 3; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.type = MTF_shape;
model_params.fiberType = 3;


% Run model
AN = modelAN(params_NT{2}, model_params); % HSR for IC input
SFIE = wrapperIC(AN.an_sout, params_NT{2}, model_params); % SFIE output

% Analysis
[rate, rate_std] = plotNaturalTimbre(params_NT{2}, SFIE.avIC, 1);

% Plot 
figure 
t = tiledlayout(3, 1, "TileSpacing","tight");
nexttile(t,1)
plot(data{1}.pitch, data{1}.rate, 'LineWidth',1.5);
hold on
plot(data{1}.pitch, rate, 'LineWidth',1.5);
title('Data and SFIE Model')
ylabel('Avg. Spike Rate (sp/s)')
xlabel('Fundamental Frequency (Hz)')


%% Model - Broad inhibition

% Model parameters
CS_param_names = {'Inhibitory strength low', 'Inhibitory strength high', 'Off-CF delay'}; % Parameters
CS_param = [0.35, 0.35, 0.003];
model_params.type = 'Lateral Model';
model_params.config_type = 'BS inhibited by off-CF BS';
model_params.lateral_CF = [CF/2 CF CF*2];
model_params.CFs = model_params.lateral_CF;
model_params.CF_range = model_params.CFs(2);
model_params.BMF = 100;
model_params.species = 1;
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 3; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.type = MTF_shape;
model_params.fiberType = 3;


% Run model
AN = modelLateralAN(params_NT{2}, model_params);
SFIE_LI = modelLateralSFIE(params_NT{2}, model_params, AN.an_sout, AN.an_sout_lo, AN.an_sout_hi,'CS_params', CS_param);

% Analysis
[rate, rate_std] = plotNaturalTimbre(params_NT{2}, SFIE_LI.avIC, 1);

% Plot 
nexttile(t,2)
plot(data{1}.pitch, data{1}.rate, 'LineWidth',1.5);
hold on
plot(data{1}.pitch, rate, 'LineWidth',1.5);
title('Data and Broad Inhibition + SFIE Model')
ylabel('Avg. Spike Rate (sp/s)')
xlabel('Fundamental Frequency (Hz)')


%% Model - Energy 

% Calculate energy prediction 
Fs = 100000;
stimulus = params_NT{2}.stim;
gamma_param.srate = Fs;
tvals = (1:length(stimulus))/Fs;
gamma_IF_reg = zeros(1,length(tvals));

%pin_gamma = zeros(size(stimulus, 1), Fs*params_NT{1}.dur+0.1*Fs);

for istim = 1:size(stimulus, 1)
	gamma_param.fc = model_params.CF; % CF
	impaired = 0; % 0 = not impaired; 1 = 'impaired'

	pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);

	pin_gamma2 = pin_gamma(istim,1:params_NT{2}.dur*Fs);
	rms_filt(istim) = sqrt(mean(pin_gamma2.^2,2));

end

% Plot 
nexttile(t,3)
plot(data{1}.pitch, data{1}.rate, 'LineWidth',1.5);
hold on
plot(data{1}.pitch, rms_filt*1000, 'LineWidth',1.5)
title('Data and Gammatone Filterbank')
ylabel('Avg. Spike Rate (sp/s)')
xlabel('Fundamental Frequency (Hz)')