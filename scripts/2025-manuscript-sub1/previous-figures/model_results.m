%% ST_model_results

%% Synthetic Timbre
clear 

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
spreadsheet_name = 'model_r2_values_ST2.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Create 

figure('position', [1796,680,7*ppi,3.2*ppi])
tiledlayout(1, 2)
fontsize = 12;
subplot_numbers = [1, 4, 7; 10, 13, 16];
labelsize = 20;

%% Simplified scatter plots 

h = gobjects(6, 1);
isSPL = sessions.SPL == 63;
for ii = 1:2

	isBS = strcmp(sessions.MTF, 'BS');
	isBE = strcmp(sessions.MTF, 'BE');
	if ii == 1
		isnotsig = sessions.p_s_e>0.05;
		x_R2 = sessions.Energy_R(isBS & ~isnotsig);
		x_R22 = sessions.Energy_R(isBE & ~isnotsig);
		x_non = sessions.Energy_R(isnotsig);
		
		y_R2 = sessions.SFIE_R(isBS & ~isnotsig);
		y_R22 = sessions.SFIE_R(isBE & ~isnotsig);
		y_non = sessions.SFIE_R(isnotsig);
	elseif ii == 2

		isnotsig = sessions.p_l_e>0.05;
		x_R2 = sessions.Energy_R(isBS & ~isnotsig);
		x_R22 = sessions.Energy_R(isBE & ~isnotsig);
		x_non = sessions.Energy_R(isnotsig);

		y_R2 = sessions.Lat_Inh_R(isBS & ~isnotsig);
		y_R22 = sessions.Lat_Inh_R(isBE & ~isnotsig);
		y_non = sessions.Lat_Inh_R(isnotsig);
	else

		isnotsig = sessions.p_l_s>0.05;
		x_R2 = sessions.SFIE_R(isBS & ~isnotsig);
		x_R22 = sessions.SFIE_R(isBE & ~isnotsig);
		x_non = sessions.SFIE_R(isnotsig);

		y_R2 = sessions.Lat_Inh_R(isBS & ~isnotsig);
		y_R22 = sessions.Lat_Inh_R(isBE & ~isnotsig);
		y_non = sessions.Lat_Inh_R(isnotsig);

	end

	h(subplot_numbers(1, ii)) = subplot(6, 3, subplot_numbers(1, ii)); % [1 4 7 10 13 16]
	hold on
	scatter(x_R2, y_R2, 20, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.7)
	scatter(x_R22, y_R22, 20, 'filled', 'MarkerEdgeColor','k', ...
		'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.7)
	%scatter(x_non, y_non, 'MarkerEdgeColor','k')
	plot([-1 1], [-1 1], 'k')
	xticklabels([])
	yticklabels([])
	set(gca, 'fontsize', fontsize)
	yline(0)
	xline(0)
	if ii == 2
		legend('BS', 'BE', '', '', 'location', 'south')
	end

	% Create distribution plot for the X-axis (horizontal)
	edges = linspace(-1, 1, 20);
	h(subplot_numbers(1, ii)+1) = subplot(6, 3, subplot_numbers(1, ii)+1);
	histogram(x_R2,edges, 'Orientation', 'vertical', 'EdgeColor', 'k');
	hold on
	histogram(x_R22,edges, 'Orientation', 'vertical', 'EdgeColor', 'k', 'FaceColor', "#D95319");
	xlim([-1 1])
	xlabel('Energy R')
	set(gca, 'fontsize', fontsize)
	xticks(-1:0.5:1)
	xline(0, 'k')
	yticks([])
	box off

	% Create distribution plot for the Y-axis (vertical) below the scatter plot
	edges = linspace(-1, 1, 20);
	h(subplot_numbers(1, ii)+2) = subplot(6, 3, subplot_numbers(1, ii)+2);
	histogram(y_R2,edges, 'Orientation', 'horizontal', 'EdgeColor', 'k');
	hold on
	histogram(y_R22,edges, 'Orientation', 'horizontal', 'EdgeColor', 'k', 'FaceColor', "#D95319");
	ylim([-1 1])
	if ii == 1
		ylabel('SFIE R')
	else
		ylabel('Broad Inhibition R')
	end
	set(gca, 'fontsize', fontsize)
	yticks(-1:0.5:1)
	yline(0, 'k')
	xticks([])
	box off
end

% Arrange plots 
all_fig_positions = ...
   [0.15,0.275,0.32,0.63;...
	0.64,0.275,0.32,0.63]; % left bottom width height

subplot_numbers = [1, 4];
for ipos = 1:2
	fig_position = all_fig_positions(ipos,:);
	nb_position = [fig_position(1),fig_position(2)-0.12,fig_position(3),0.1];
	wb_position = [fig_position(1)-0.05,fig_position(2),0.04,fig_position(4)];
	set(h(subplot_numbers(ipos)), 'Position', fig_position)
	set(h(subplot_numbers(ipos)+1), 'Position', nb_position)
	set(h(subplot_numbers(ipos)+2), 'Position', wb_position)
end

% Annotate
labels = {'A', 'B',};
labelleft= repmat([0.01, 0.51], 1, 2);
labelbottom = [repmat(0.95,1, 2) repmat(0.48, 1, 2)];
for ii = 1:2
	annotation('textbox',[labelleft(ii) labelbottom(ii) 0.071 0.058],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% RMSE Plots 
% 
% figure('position', [1796,680,7*ppi,3.2*ppi])
% fontsize = 12;
% subplot_numbers = [1, 4, 7; 10, 13, 16];
% labelsize = 20;
% 
% h = gobjects(6, 1);
% isSPL = sessions.SPL == 63;
% for ii = 1:2
% 
% 	isBS = strcmp(sessions.MTF, 'BS');
% 	isBE = strcmp(sessions.MTF, 'BE');
% 	if ii == 1
% 		x_R2 = sessions.Energy_rmse(isBS);
% 		x_R22 = sessions.Energy_rmse(isBE);
% 
% 		y_R2 = sessions.SFIE_rmse(isBS);
% 		y_R22 = sessions.SFIE_rmse(isBE);
% 	elseif ii == 2
% 		x_R2 = sessions.Energy_rmse(isBS);
% 		x_R22 = sessions.Energy_rmse(isBE);
% 
% 		y_R2 = sessions.Lat_Inh_rmse(isBS);
% 		y_R22 = sessions.Lat_Inh_rmse(isBE);
% 	else
% 		x_R2 = sessions.SFIE_rmse(isBS);
% 		x_R22 = sessions.SFIE_rmse(isBE);
% 
% 		y_R2 = sessions.Lat_Inh_rmse(isBS);
% 		y_R22 = sessions.Lat_Inh_rmse(isBE);
% 	end
% 
% 	h(subplot_numbers(1, ii)) = subplot(6, 3, subplot_numbers(1, ii)); % [1 4 7 10 13 16]
% 	hold on
% 	scatter(x_R2, y_R2, 20, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.7)
% 	scatter(x_R22, y_R22, 20, 'filled', 'MarkerEdgeColor','k', ...
% 		'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.7)
% 	%scatter(x_non, y_non, 'MarkerEdgeColor','k')
% 	plot([-1 90], [-1 90], 'k')
% 	xticks(0:20:80)
% 	yticks(0:20:80)
% 	xticklabels([])
% 	yticklabels([])
% 	set(gca, 'fontsize', fontsize)
% 	yline(0)
% 	xline(0)
% 	ylim([0 90])
% 	xlim([0 90])
% 	if ii == 1
% 		legend('BS', 'BE', '', '', 'location', 'south')
% 	end
% 	grid on
% 
% 	% Create distribution plot for the X-axis (horizontal)
% 	edges = linspace(0, 90, 20);
% 	h(subplot_numbers(1, ii)+1) = subplot(6, 3, subplot_numbers(1, ii)+1);
% 	histogram(x_R2,edges, 'Orientation', 'vertical', 'EdgeColor', 'k');
% 	hold on
% 	histogram(x_R22,edges, 'Orientation', 'vertical', 'EdgeColor', 'k', 'FaceColor', "#D95319");
% 	xlim([0 90])
% 	xlabel('Energy RMSE_n_o_r_m')
% 	set(gca, 'fontsize', fontsize)
% 	xticks(0:20:80)
% 	yticks([])
% 	grid on
% 	box off
% 
% 	% Create distribution plot for the Y-axis (vertical) below the scatter plot
% 	edges = linspace(0, 90, 20);
% 	h(subplot_numbers(1, ii)+2) = subplot(6, 3, subplot_numbers(1, ii)+2);
% 	histogram(y_R2,edges, 'Orientation', 'horizontal', 'EdgeColor', 'k');
% 	hold on
% 	histogram(y_R22,edges, 'Orientation', 'horizontal', 'EdgeColor', 'k', 'FaceColor', "#D95319");
% 	ylim([0 90])
% 	if ii == 1
% 		ylabel('SFIE RMSE_n_o_r_m')
% 	else
% 		ylabel('Broad Inh. RMSE_n_o_r_m')
% 	end
% 	set(gca, 'fontsize', fontsize)
% 	yticks(0:20:80)
% 	xticks([])
% 	grid on
% 	box off
% end
% 
% % Arrange plots 
% all_fig_positions = ...
%    [0.15,0.29,0.32,0.63;...
% 	0.64,0.29,0.32,0.63]; % left bottom width height
% 
% subplot_numbers = [1, 4];
% for ipos = 1:2
% 	fig_position = all_fig_positions(ipos,:);
% 	nb_position = [fig_position(1),fig_position(2)-0.12,fig_position(3),0.1];
% 	wb_position = [fig_position(1)-0.05,fig_position(2),0.04,fig_position(4)];
% 	set(h(subplot_numbers(ipos)), 'Position', fig_position)
% 	set(h(subplot_numbers(ipos)+1), 'Position', nb_position)
% 	set(h(subplot_numbers(ipos)+2), 'Position', wb_position)
% end
% 
% % Annotate
% labels = {'C', 'D',};
% labelleft= repmat([0.01, 0.51], 1, 2);
% labelbottom = [repmat(0.95,1, 2) repmat(0.48, 1, 2)];
% for ii = 1:2
% 	annotation('textbox',[labelleft(ii) labelbottom(ii) 0.071 0.058],...
% 		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
% 		'EdgeColor','none');
% end