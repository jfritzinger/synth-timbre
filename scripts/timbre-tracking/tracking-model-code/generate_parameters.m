function [Fc_high, Fc_low, F0_high, F0_low] = generate_parameters(task, nom_Fc, nom_F0, delta, rove, DL_var)
% Returns the geometrically centered spectral centroid and F0 for
% experiments 

%% Determine difference limen for the task

if task == 1 % timbre
    DL = 0.8; % pitch DL, musicians Allen & Oxenham 2014
else % pitch
    DL = 4; % timbre DL, musicians Allen & Oxenham 2014
end

%% Calculate F0 and Fc values 

switch task
    case 1 % timbre experiment 
        
        % Set timbre values
        if rove == 2  % calculate rove and geometrically center around nominal
            delta_Fc = 1 + delta;
            nom_Fc_roved = nom_Fc + (rand*nom_Fc*0.2 - nom_Fc*0.1);  % roved nominal Fc
            Fc_low = sqrt((nom_Fc_roved^2)/delta_Fc);   % calculate geomean low Fc
            Fc_high = Fc_low*delta_Fc;                  % calculate geomean high Fc
        else % no rove
            delta_Fc = 1 + delta;
            Fc_low = sqrt((nom_Fc^2)/delta_Fc);
            Fc_high = Fc_low*delta_Fc;
        end
        
        % Set pitch values
        if DL_var == -1 % no variation or rove
            F0_high = nom_F0;
            F0_low = nom_F0;
        else % calculate rove and geometrically center around nominal
            nom_F0_roved = nom_F0 + (rand*nom_F0*0.2 - nom_F0*0.1); % roved nominal value
            DL_change = nom_F0*DL*DL_var;
            delta_F0 = 1 + DL_change/nom_F0;
            F0_low = sqrt((nom_F0_roved^2)/delta_F0);
            F0_high = F0_low*delta_F0;
        end
        
    case 2 % pitch experiment 
        
        % Set pitch values 
        if rove == 2  % calculate rove and geometrically center around nominal
            delta_F0 = 1 + delta;
            nom_F0_roved = nom_F0 + (rand*nom_F0*0.2 - nom_F0*0.1);  % roved nominal F0
            F0_low = sqrt((nom_F0_roved^2)/delta_F0);   % calculate geomean low F0
            F0_high = F0_low*delta_F0;                  % calculate geomean high F0
        else % no rove
            delta_F0 = 1 + delta;
            F0_low = sqrt((nom_F0^2)/delta_F0);
            F0_high = F0_low*delta_F0;
        end
        
        % Set timbre values
        if DL_var == -1 % no variation or rove
            Fc_high = nom_Fc;
            Fc_low = nom_Fc;
        else % calculate rove and geometrically center around nominal
            nom_Fc_roved = nom_Fc + (rand*nom_Fc*0.2 - nom_Fc*0.1); % roved nominal value
            DL_change = nom_Fc*DL*DL_var;
            delta_Fc = 1 + DL_change/nom_Fc;
            Fc_low = sqrt((nom_Fc_roved^2)/delta_Fc);
            Fc_high = Fc_low*delta_Fc;
        end

end

