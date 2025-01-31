%% Synthetic Timbre
clear 

% Load in spreadsheet
[base, datapath, ~, ~] = getPaths();
spreadsheet_name = 'model_r2_values_ST.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Simplified scatter plots 

figure('position', [1796,680,1292,321])
tiledlayout(1, 4)
spls = [43, 63, 73, 83];
ispl = 2;
isSPL = sessions.SPL == spls(ispl);

isMTF = strcmp(sessions.MTF, 'BS');
BS_ind = find(isMTF & isSPL);
SFIE_R2 = sessions.SFIE_R(BS_ind);
energy_R2 = sessions.Energy_R(BS_ind);
lat_inh_R2 = sessions.Lat_Inh_R(BS_ind);
pop_R2 = sessions.SFIE_Pop_R(BS_ind);

isMTF = strcmp(sessions.MTF, 'BE');
BE_ind = find(isMTF & isSPL);
SFIE_R22 = sessions.SFIE_R(BE_ind);
energy_R22 = sessions.Energy_R(BE_ind);
lat_inh_R22 = sessions.Lat_Inh_R(BE_ind);
pop_R22 = sessions.SFIE_Pop_R(BE_ind);

nexttile
hold on
scatter(energy_R2, SFIE_R2, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.5)
scatter(energy_R22, SFIE_R22, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.5)
plot([-1 1], [-1 1], 'k')
xlabel('Energy R')
ylabel('SFIE R')
title('SFIE / Energy')
set(gca, 'fontsize', 16)
yline(0)
xline(0)

nexttile
hold on
scatter(energy_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.5)
scatter(energy_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k',...
	'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.5)
plot([-1 1], [-1 1], 'k')
xlabel('Energy R')
ylabel('Lateral Inh R')
set(gca, 'fontsize', 16)
yline(0)
xline(0)
title('Lat Inh / Energy')

nexttile
hold on
scatter(SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.5)
scatter(SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.5)
plot([-1 1], [-1 1], 'k')
xlabel('SFIE R')
ylabel('Lateral Inh R')
set(gca, 'fontsize', 16)
yline(0)
xline(0)
legend('BS', 'BE', '', '')
title('Lat Inh / SFIE')

nexttile
hold on
scatter(SFIE_R2, pop_R2, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.5)
scatter(SFIE_R22, pop_R22, 'filled', 'MarkerEdgeColor','k',...
	'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.5)
plot([-1 1], [-1 1], 'k')
xlabel('SFIE R')
ylabel('SFIE population R')
set(gca, 'fontsize', 16)
yline(0)
xline(0)
title('SFIE / SFIE Population')

%% Histograms 

%Plot histograms 
figure('position', [47,324,523,564])
tiledlayout(3, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 2

	MTF_target = 'BS';
	isMTF = strcmp(sessions.MTF, MTF_target);
	isSPL = sessions.SPL == spls(ispl);
	indices = find(isMTF & isSPL);

	nexttile
	hold on
	CFs = sessions.CF(indices);
	SFIE_R2 = sessions.SFIE_R(indices);
	energy_R2 = sessions.Energy_R(indices);
	lat_inh_R2 = sessions.Lat_Inh_R(indices);
	edges = linspace(-1, 1, 20);
	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)
	title(['BS SFIE, ' num2str(spls(ispl)) ' dB SPL'])

	nexttile
	histogram(energy_R2, edges, 'FaceColor','#D55E00', 'FaceAlpha',0.5)
	title('BS, energy')

	nexttile
	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
	title('BS, Lateral Inhibition')

	MTF_target = 'BE';
	isMTF = strcmp(sessions.MTF, MTF_target);
	isSPL = sessions.SPL == spls(ispl);
	indices = find(isMTF & isSPL);

	nexttile
	hold on
	CFs = sessions.CF(indices);
	SFIE_R2 = sessions.SFIE_R(indices);
	energy_R2 = sessions.Energy_R(indices);
	lat_inh_R2 = sessions.Lat_Inh_R(indices);
	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)
	title(['BE SFIE, ' num2str(spls(ispl)) ' dB SPL'])
	
	nexttile
	histogram(energy_R2, edges, 'FaceColor','#D55E00', 'FaceAlpha',0.5)
	title('BE, energy')

	nexttile
	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
	title('BE, Lateral Inhibition')
end