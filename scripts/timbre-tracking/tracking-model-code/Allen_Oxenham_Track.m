%function Allen_Oxenham_Track()
%% Allen and Oxenham Track
% Johanna Fritzinger
% 7/9/2020
clear all
close all

%% Naming Conventions

con_name = {'Timbre tracking' 'Pitch tracking'};
rove_name = {'no rove' 'rove'};
DL_var_name = {'0DL, no rove' '0DL' '2DL' '5DL' '10DL' '25DL' '50DL' '100DL'};

%% Starting Parameters

% General conditions
DL_var_val = [-1, 0, 0.02, 0.05, 0.1, 0.25, 0.5, 1]; % variation in non-target dim, -1 is 0DL no rove
num_reps = 1;       % number of reps
condition = 1;      % 1 = timbre, 2 = pitch
rove = 2;           % 1 = no rove, 2 = rove
DL_variation = 1:8;   % variation in non-target dimension
plt = 0;           % 0 = no plot, 1 = plot trials 

% Timbre condition
delta_Fc = 0.25;    % Hz, 25% CF
Fc_nom = 1200;      % Hz, starting value

% Pitch condition
delta_F0 = 0.25;    % Hz, 10% F0
F0_nom = 200;       % Hz, starting value

% Spectral centroid parameters
dur = 0.5;          % sec, stimulus duration
rampdur = 0.02;     % sec, ramp duration
isi = 0.3;          % sec, interstimulus interval
G = 24;             % dB/octave, slope of triangle edges
stimdB = 70;        % dB, overall level
Fs = 100e3;         % model sampling rate

%% Tracking
% Based on the methods from Allen & Oxenham 2014

for n_reps = num_reps
    for con = condition
        for r = rove
            for DL_var = DL_variation
                
                % Initalize tracking parameters
                nreversals = 0; % counts number of reversals
                ncorrect = 0; % 0-2 based on how many trials correct in a row
                reps = 1; % counts number of trials
                direction = [1 1];
                reversal_tracker = 0; % records what trials reversals occured
                delta_factor = 2;
                delta_tracker = zeros(100, 1);
                
                % Set values for pitch or timbre condition
                if con == 1
                    delta = delta_Fc;
                else
                    delta = delta_F0;
                end
                
                while nreversals < 8
                    delta_tracker(reps) = delta;
                    
                    % Create Fc and F0 parameters
                    [Fc_high, Fc_low, F0_high, F0_low] = generate_parameters(con, Fc_nom, F0_nom, delta, r, DL_var_val(DL_var));
                    
                    % Create stimulus
                    F0 = [F0_low, F0_high];
                    Fc = [Fc_low, Fc_high];
                    F0_order = randperm(length(F0)); % randomize order of F0
                    Fc_order = randperm(length(Fc)); % randomize order of Fc
                    stim(1,:) = generate_spectral_centroid(dur, rampdur, F0(F0_order(1)), Fc(Fc_order(1)), G, stimdB, Fs);
                    stim(2,:) = generate_spectral_centroid(dur, rampdur, F0(F0_order(2)), Fc(Fc_order(2)), G, stimdB, Fs);
                    
                    % Call model
                    stimulus = [stim(1,:), zeros(1, Fs*isi), stim(2, :)];
                    %DV_human(reps) = play_experiment(stimulus, Fs, Fc_order);
                    [DV_model(reps), est_Fc(reps,:)] = spec_centroid_model(condition, stim, Fc, F0, Fc_order, F0_order, plt);
                    
                    % Tracking
                    if DV_model(reps) == 1 % correct trial
                        ncorrect = ncorrect+1;
                    else % incorrect trial
                        ncorrect = 0;
                        direction(1) = direction(2);
                        direction(2) = 0;
                        delta = delta*delta_factor;
                    end
                    if ncorrect == 2 % decrease delta after 2 correct trials
                        delta = delta/delta_factor;
                        ncorrect = 0;
                        direction(1) = direction(2);
                        direction(2) = 1;
                    elseif ncorrect == 1
                        direction(1) = direction(2);
                    end
                    if direction(2) ~= direction(1)
                        nreversals = nreversals +1;
                        reversal_tracker(reps) = 1;
                    else
                        reversal_tracker(reps) = 0;
                    end
                    if nreversals == 2
                        delta_factor = 1.26;
                    elseif nreversals == 4
                        delta_factor = 1.12;
                    end
                    reps = reps+1;
                end
                delta_tracker = nonzeros(delta_tracker);
                percentage = delta_tracker.*100;
                reversal_tracker(reversal_tracker == 0) = NaN;
                DV_model(DV_model == 0) = NaN;
                DV_model(DV_model == 1) = 0;
                
                % Plot tracking
                figure;
                %subplot(2,1,1);
                plot(percentage);
                hold on
                plot(reversal_tracker.*percentage', 'X', 'MarkerSize', 7)
                %plot(DV_model, '-r.', 'MarkerSize', 13)
                title([con_name{con} ' with ' rove_name{r} ', ' DL_var_name{DL_var} ' in non-target dimension'])
                ylabel('Threshold Percentage')
                xlabel('Trials')
                %legend('Human track', 'Human reversals', 'Model errors');
                
                % Plot model
                %                 subplot(2,1,2);
                %                 plot(est_Fc, '.', 'MarkerSize', 13);
                %                 title('Estimated Model Spectral Centroids, Smoothes, Subtracted CoM DV')
                
                % Psychometric function
                reversal_tracker(isnan(reversal_tracker)) = 0;
                start_ind = find(reversal_tracker);
                thresh(DL_var) = mean(percentage(start_ind(3):start_ind(end)));
            end
        end
    end
end
%end

function DV = play_experiment(stim, Fs, order)

% start experiment
disp('Listen for the tone with brighter timbre. Press 1 for the first tone, press 2 for the second tone.');
sound(stim, Fs);
x = input('Answer: ');

% evaluate decision
if order(1) == 1 % low played first
    if x == 2
        DV = 1;
        disp('Correct');
    elseif x == 1
        DV = 0;
        disp('Incorrect');
    end
else % high played first
    if x == 2
        DV = 0;
        disp('Incorrect');
    elseif x == 1
        DV = 1;
        disp('Correct');
    end
end
end