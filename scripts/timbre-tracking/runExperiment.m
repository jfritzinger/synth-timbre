function DV = runExperiment(stim, Fs, order)
% Plays two tones and makes the user choose the brighter or higher tone
% J. Fritzinger, updated 1/6/2021

% Start experiment
disp('Listen for the tone with brighter timbre. Press 1 for the first tone, press 2 for the second tone.');
sound(stim, Fs);
x = input('Answer: ');

% Evaluate decision
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