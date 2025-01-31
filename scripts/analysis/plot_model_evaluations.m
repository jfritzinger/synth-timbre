%% plot_model_evaluations.m
%
% Script that reads in table from create_r2_table.m with all variance
% explained metrics for each model. Used to create plots evaluating which
% models most correctly predict the physiological responses, in general and
% over different CF ranges. 
%
%
% Author: J. Fritzinger
% Created: 2022-09-13; Last revision: 2024-09-26
%
% -------------------------------------------------------------------------
clear

%% Synthetic Timbre

% Load in spreadsheet
[base, datapath, ~, ~] = getPaths();
spreadsheet_name = 'model_r2_values_ST.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Plot variance explained by the energy model (ST)

% Find sessions for target MTF type
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);

figure('position', [519,592,1100,240])
tiledlayout(2, 4, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 1:4
	isSPL = sessions.SPL == spls(ispl);
	indices = find(isMTF & isSPL);

	nexttile
	hold on
	CFs = sessions.CF(indices);
	SFIE_R2 = sessions.SFIE_R2(indices);
	energy_R2 = sessions.Energy_R2(indices);
	scatter(CFs, SFIE_R2, 20, 'filled', 'MarkerFaceColor','#009E73')
	set(gca, 'XScale', 'log')

	nexttile
	hold on
	scatter(CFs, energy_R2, 20, 'filled', 'MarkerFaceColor','#D55E00')
	set(gca, 'XScale', 'log')
end

%%

%Plot histograms 
figure('position', [47,487,523,401])
tiledlayout(6, 4, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 1:4

	MTF_target = 'BS';
	isMTF = strcmp(sessions.MTF, MTF_target);
	isSPL = sessions.SPL == spls(ispl);
	indices = find(isMTF & isSPL);

	nexttile
	hold on
	CFs = sessions.CF(indices);
	SFIE_R2 = sessions.SFIE_R2(indices);
	energy_R2 = sessions.Energy_R2(indices);
	lat_inh_R2 = sessions.Lat_Inh_R2(indices);
	edges = linspace(0, 1, 20);
	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)

	title(['BS SFIE, ' num2str(spls(ispl)) ' dB SPL'])
	nexttile
	hold on
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
	SFIE_R2 = sessions.SFIE_R2(indices);
	energy_R2 = sessions.Energy_R2(indices);
	lat_inh_R2 = sessions.Lat_Inh_R2(indices);
	edges = linspace(0, 1, 20);
	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)
	title(['BE SFIE, ' num2str(spls(ispl)) ' dB SPL'])
	
	nexttile
	hold on
	histogram(energy_R2, edges, 'FaceColor','#D55E00', 'FaceAlpha',0.5)
	title('BE, energy')
	%title(['BE, ' num2str(spls(ispl)) ' dB SPL'])

	nexttile
	hold on
	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
	title('BE, Lateral Inhibition')
end

%% Scatter plots 

% energy vs SFIE
figure('position', [78,337,939,616])
tiledlayout(3, 4, 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 1:4
	isSPL = sessions.SPL == spls(ispl);

	isMTF = strcmp(sessions.MTF, 'BS');
	BS_ind = find(isMTF & isSPL);
	SFIE_R2 = sessions.SFIE_R(BS_ind);
	energy_R2 = sessions.Energy_R(BS_ind);
	lat_inh_R2 = sessions.Lat_Inh_R(BS_ind);

	isMTF = strcmp(sessions.MTF, 'BE');
	BE_ind = find(isMTF & isSPL);
	SFIE_R22 = sessions.SFIE_R(BE_ind);
	energy_R22 = sessions.Energy_R(BE_ind);
	lat_inh_R22 = sessions.Lat_Inh_R(BE_ind);

	nexttile
	hold on
	scatter(energy_R2, SFIE_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, SFIE_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, SFIE_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(energy_R22, SFIE_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('Energy R')
	ylabel('SFIE R')
	title(sprintf('%d dB SPL', spls(ispl)))

	nexttile
	hold on
	scatter(energy_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, lat_inh_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(energy_R22, lat_inh_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('Energy R')
	ylabel('Lateral Inh R')

	nexttile
	hold on
	scatter(SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(SFIE_R2, lat_inh_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(SFIE_R22, lat_inh_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('SFIE R')
	ylabel('Lateral Inh R')
	legend('BS', 'BE')
end

%% 3D Scatter Plot 
% 
% figure('position', [78,337,939,616])
% tiledlayout(1, 4, 'TileIndexing','columnmajor')
% spls = [43, 63, 73, 83];
% for ispl = 1:4
% 	isSPL = sessions.SPL == spls(ispl);
% 
% 	isMTF = strcmp(sessions.MTF, 'BS');
% 	BS_ind = find(isMTF & isSPL);
% 	SFIE_R2 = sessions.SFIE_R(BS_ind);
% 	energy_R2 = sessions.Energy_R(BS_ind);
% 	lat_inh_R2 = sessions.Lat_Inh_R(BS_ind);
% 
% 	isMTF = strcmp(sessions.MTF, 'BE');
% 	BE_ind = find(isMTF & isSPL);
% 	SFIE_R22 = sessions.SFIE_R(BE_ind);
% 	energy_R22 = sessions.Energy_R(BE_ind);
% 	lat_inh_R22 = sessions.Lat_Inh_R(BE_ind);
% 
% 	nexttile
% 	hold on
% 	scatter3(energy_R2, SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
% 	scatter3(energy_R22, SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
% 	xlabel('Energy R^2')
% 	ylabel('SFIE R^2')
% 	zlabel('Lat Inh R^2')
% 	title(sprintf('%d dB SPL', spls(ispl)))
% 	grid on
% 
% end

%% Simplified scatter plots 

figure('position', [78,337,939,616])
tiledlayout(1, 3, 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 2
	isSPL = sessions.SPL == spls(ispl);

	isMTF = strcmp(sessions.MTF, 'BS');
	BS_ind = find(isMTF & isSPL);
	SFIE_R2 = sessions.SFIE_R(BS_ind);
	energy_R2 = sessions.Energy_R(BS_ind);
	lat_inh_R2 = sessions.Lat_Inh_R(BS_ind);

	isMTF = strcmp(sessions.MTF, 'BE');
	BE_ind = find(isMTF & isSPL);
	SFIE_R22 = sessions.SFIE_R(BE_ind);
	energy_R22 = sessions.Energy_R(BE_ind);
	lat_inh_R22 = sessions.Lat_Inh_R(BE_ind);

	nexttile
	hold on
	scatter(energy_R2, SFIE_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, SFIE_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, SFIE_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(energy_R22, SFIE_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('Energy R')
	ylabel('SFIE R')
	title(sprintf('%d dB SPL', spls(ispl)))
	set(gca, 'fontsize', 16)

	nexttile
	hold on
	scatter(energy_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, lat_inh_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(energy_R22, lat_inh_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('Energy R')
	ylabel('Lateral Inh R')
	set(gca, 'fontsize', 16)

	nexttile
	hold on
	scatter(SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(SFIE_R2, lat_inh_R2);
	%plot(x, y, 'b')
	[x, y] = findRegressionLine(SFIE_R22, lat_inh_R22);
	plot([-1 1], [-1 1], 'k')
	%plot(x, y, 'r')
	xlabel('SFIE R')
	ylabel('Lateral Inh R')
	legend('BS', 'BE')
	set(gca, 'fontsize', 16)
end

%% Functions 

function [x, y] = findRegressionLine(datax, datay)
	x = -1:0.05:1;


	mdl = fitlm(datax, datay);
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	y = p(1)*x+p(2);

end

