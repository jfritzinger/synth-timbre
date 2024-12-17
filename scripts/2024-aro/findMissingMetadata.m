%% Checking stimulus_sessions datasheet for missing metadata to fill out

%% Reads in stimuli_sessions

if ismac
    mount('nsc-lcarney-g1')
    %filename = '/Volumes/Rabbit_data/session_table.xlsx';
    filename = '/Volumes/Rabbit_data/session_table.xlsx';
else
    filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
end

sessions = readtable(filename, 'PreserveVariableNames',true);

%% Check sessions for a certain stimulus 

rabbit = [24, 25, 26, 27, 28, 29, 30];
rating = {'Good', 'Excellent'};

%stimulus2 = 'SPEC_slide_WB_noise';
stimulus = 'SPEC_slide_Spectral_Centroid';
%stimulus = 'click_ITD';

%% Gets sessions of interest

interest_rabbit = ismember(sessions.Rabbit,rabbit);
interest_rating = ismember(sessions.Rating,rating);

names = sessions.Properties.VariableNames;
ind = find(strcmp(names,stimulus));
stim_columns = table2array(sessions(:,ind));
stim_columns(stim_columns > 1) = 1; % ignores Paul's system of using 1 and 2
interest_stim = sum(stim_columns,2);

% ind = find(strcmp(names,stimulus2));
% stim_columns = table2array(sessions(:,ind));
% stim_columns(stim_columns > 1) = 1; % ignores Paul's system of using 1 and 2
% interest_stim2 = sum(stim_columns,2);

%% Checks sessions for missing data 

%empty_depth = isnan(sessions.Tetrode_Depth);
empty_CF = isnan(sessions.CF);
%empty_MTF = isnan(sessions.MTF);

%interest = interest_rabbit & interest_rating & empty_CF & (interest_stim |interest_stim2);
interest = interest_rabbit & interest_rating & empty_CF & interest_stim;

num_interest = sum(interest);
index = find(interest);

session_interest = {};
for ii = 1:num_interest
    session_interest{ii} = sprintf('R0%dS%d TT%d N%d', sessions.Rabbit(index(ii)), ...
        sessions.Session(index(ii)), sessions.Tetrode(index(ii)), ...
        sessions.Neuron(index(ii)));
end

%% Create and save table

addData = table('Size',[500 3], 'VariableTypes', ["string", "double", "double"],...
    'VariableNames',["Session", "CF", "Depth"]);
for ii = 1:num_interest
    addData.Session(ii) = session_interest(ii);
    addData.CF(ii) = sessions.CF(index(ii));
    %addData.Depth(ii) = sessions.Tetrode_Depth(index(ii));
end

writetable(addData,'findMetadataToAdd.xls')

