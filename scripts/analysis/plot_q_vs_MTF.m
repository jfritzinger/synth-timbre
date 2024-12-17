%% plot_q_vs_MTFslope_example.m
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking2.xlsx"));


%% Analysis

figure('position', [519,299,1150,539])
tiledlayout(1, 4, 'TileSpacing','tight', 'Padding','compact')

isbin = tables.binmode == 2;
is200 = tables.F0 == 200;
spls = [43, 63, 73, 83];
for ispl = 1:4

	% Get data
	islevel = tables.SPL == spls(ispl);
	index = islevel & isbin & is200;

	% Data
	CFs = tables.CF(index);
	Qs = tables.Q(index);
	MTFs = tables.MTF(index);
	MTF_str = tables.MTF_str(index);

	% Plot
	nexttile
	gscatter(MTF_str, Qs, MTFs, 'filled')
	hold on
	box on

end

%% 

ispeak = strcmp(tables.Type, 'Peak');
isdip = strcmp(tables.Type, 'Dip');
isflat = strcmp(tables.Type, 'Flat');
signed_Q = tables.Q;
signed_Q(isdip) = signed_Q(isdip) * -1;
signed_Q(isflat) = 0;

isBE = strcmp(tables.MTF, 'BE');
isBS = strcmp(tables.MTF, 'BS');
%ishybrid = contains(tables.MTF, 'H');
%isflat = strcmp(tables.MTF, 'F');
signed_MTF = tables.MTF_str;
% signed_MTF(isBE) = signed_MTF(isBE) * -1;

isbin = tables.binmode == 2;
is200 = tables.F0 == 200;
spls = [43, 63, 73, 83];
figure('position', [73,548,1596,290])
tiledlayout(1, 4, 'TileSpacing','tight', 'Padding','compact')
for ispl = 1:4

	% Get data
	islevel = tables.SPL == spls(ispl);
	index = islevel & isbin & is200;
	isBE = strcmp(tables.MTF(index), 'BE');
	isBS = strcmp(tables.MTF(index), 'BS');
	

	% Data
	CFs = tables.CF(index);
	Qs = signed_Q(index);
	MTFs = tables.MTF(index);
	MTF_str = signed_MTF(index);

	% Fit lines 

	% Plot
	nexttile
	hold on
	scatter(MTF_str(isBS), Qs(isBS), 'filled', 'MarkerEdgeColor','k')
	scatter(MTF_str(isBE), Qs(isBE), 'filled', 'MarkerEdgeColor','k')
	%gscatter(MTF_str, Qs, MTFs, 'filled')
	xline(0)
	yline(0)
	box on
	ylabel('Q')
	xlabel('MTF % Change')
	title(sprintf('%d dB SPL', spls(ispl)))
	legend('BS', 'BE')

end


%% Plot
% MTF = categorical(MTF);
% figure
% hold on
% xline(0)
% yline(0)
% gscatter(slope, q, MTF, 'filled')
% % xlim([-7 7])
% % ylim([-7 7])
% 
% xlabel('MTF max slope')
% ylabel('Q_2_5')
% title('Quantification of peak sharpness vs MTF slope')