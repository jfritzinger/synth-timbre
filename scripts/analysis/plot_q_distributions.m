%% plot_q_distributions.m
%
% Script that loads in spreadsheet of q values and plots them based on
%
%
% Author: J. Fritzinger
% Created: 2022-10-21; Last revision: 2024-10-21
%
% -------------------------------------------------------------------------
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking.xlsx"));

%% Plot distribution of Q values

figure

q_vals = tables{1:end,14};
q_vals(isnan(q_vals)) = [];

edges = linspace(0, 30, 30);
histogram(q_vals, edges)

lims = tables{1:end,12};
num_data = sum(~isnan(lims));
num_lims = sum(lims, 'omitnan');

title(sprintf('%s, %d/%d (%0.0f percent) hit limits', 'Q', num_lims, num_data, num_lims/num_data*100))


%% Plot contra vs binaural Q values 

figure
isbin = tables.binmode == 2;
iscontra = tables.binmode == 1;

% bin
q_vals = tables{isbin,14};
q_vals(isnan(q_vals)) = [];
edges = linspace(0, 50, 50);
histogram(q_vals, edges)
hold on

% contra
q_vals = tables{iscontra,14};
q_vals(isnan(q_vals)) = [];
edges = linspace(0, 50, 50);
histogram(q_vals, edges)

title('Q')
legend('bin', 'contra')

%% Plot different sound level Q values 

figure
is43 = tables.SPL == 43;
is63 = tables.SPL == 63;
is83 = tables.SPL == 83;

% 43 dB SPL 
q_vals = tables{is43,14};
q_vals(isnan(q_vals)) = [];
edges = linspace(0, 50, 50);
histogram(q_vals, edges)
hold on

% 43 dB SPL 
q_vals = tables{is63,14};
q_vals(isnan(q_vals)) = [];
edges = linspace(0, 50, 50);
histogram(q_vals, edges)
hold on

% 43 dB SPL 
q_vals = tables{is83,14};
q_vals(isnan(q_vals)) = [];
edges = linspace(0, 50, 50);
histogram(q_vals, edges)
hold on

title('Q')
legend('43 dB SPL', '63 dB SPL', '83 dB SPL')

