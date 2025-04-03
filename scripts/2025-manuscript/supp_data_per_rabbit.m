%% Fig0_DistributionOfData
% J. Fritzinger, updated 12/15/23
%
% This script loads in the putative neurons spreadsheet and plots the MTF
% distribution, BMF distribution, WMF distribution, hybrid BMF/WMF
% distribution, and CF distribution for all neurons combined and split into
% the four rabbits (R24, R25, R27, R29)
clear

%% Load in spreadsheet

[base, ~, savepath, ppi] = getPaths();
datapath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_units = size(sessions, 1);
color = {"#0072BD", "#D95319", '#EDB120', '#7E2F8E'};

%% Get number of putative neurons from each rabbit

figure('Position',[247,154,1294,1002])
tiledlayout(6, 5, "TileSpacing","tight", 'Padding','compact')
font_size = 14;

num_putative_units = size(sessions, 1);
putative_units = sessions.Putative_Units;
rabbit_units = cellfun(@(p) str2double(p(2:3)), putative_units(1:334), 'UniformOutput', false);
rabbit_putative_units = [rabbit_units{:}];
rabbit = [24, 25, 27, 29];
rabbit_ages = {'44-54 months', '44-61 months', '6-29 months', '~18-33 months'};

number_units = zeros(4, 1);
for irab = 1:4
	rab = rabbit(irab);
	num_units = sum(rabbit_putative_units==rab);
	number_units(irab) = num_units;
end

nexttile([1 5])
rabbits = ["R024";"R025";"R027";"R029"];
if ismac
	piechart(number_units, rabbits)
else
	pie(number_units, rabbits)
end
title(['Putative Neurons per Rabbit, n = ' num2str(num_putative_units)])
set(gca, 'fontsize', font_size)

%% Get MTFs for each putative neuron
% WBTIN Diotic

MTFs = sessions.MTF;
num_sesh = length(MTFs);
MTF_type = zeros(num_sesh,1);
for isesh = 1:num_sesh
	MTF_shape = MTFs{isesh};
	if contains(MTF_shape, 'H')
		MTF_type(isesh) = 3;
	elseif strcmp(MTF_shape, 'BE')
		MTF_type(isesh) = 1;
	elseif strcmp(MTF_shape, 'BS')
		MTF_type(isesh) = 2;
	else % Flat
		MTF_type(isesh) = 4;
	end
end

MTF_names = categorical({'BE','BS','Hybrid','Flat'});
MTF_names = reordercats(MTF_names,{'BE','BS','Hybrid','Flat'});
num_BE = sum(MTF_type==1);
num_BS = sum(MTF_type==2);
num_H = sum(MTF_type==3);
num_n = sum(MTF_type==4);
num_types = [num_BE num_BS num_H num_n];

nexttile()
bar(MTF_names,num_types, 'black');
ylabel('# Neurons')
title(sprintf('MTF Types, n=%d', sum(num_types)))
set(gca, 'FontSize', font_size)
grid on
box on

% MTF per rabbits
putative_units = sessions.Putative_Units;
rabbit_units = cellfun(@(p) str2double(p(2:3)), putative_units(1:334), 'UniformOutput', false);
rabbit_putative_units = [rabbit_units{:}];
rabbit = [24, 25, 27, 29];
for irab = 1:4
	clear MTF_type

	rab = rabbit(irab);
	has_rab = rabbit_putative_units==rab;
	MTFs = sessions.MTF(has_rab);
	num_sesh = length(MTFs);
	MTF_type = zeros(num_sesh,1);
	for isesh = 1:num_sesh
		MTF_shape = MTFs{isesh};
		if contains(MTF_shape, 'H')
			MTF_type(isesh) = 3;
		elseif strcmp(MTF_shape, 'BE')
			MTF_type(isesh) = 1;
		elseif strcmp(MTF_shape, 'BS')
			MTF_type(isesh) = 2;
		else % Flat
			MTF_type(isesh) = 4;
		end
	end

	MTF_names = categorical({'BE','BS','Hybrid','Flat'});
	MTF_names = reordercats(MTF_names,{'BE','BS','Hybrid','Flat'});
	num_BE = sum(MTF_type==1);
	num_BS = sum(MTF_type==2);
	num_H = sum(MTF_type==3);
	num_n = sum(MTF_type==4);
	num_types = [num_BE num_BS num_H num_n];

	nexttile
	bar(MTF_names,num_types, 'FaceColor', color{irab}, 'EdgeColor',color{irab});
	title(sprintf('R0%d, %s, n=%d', rabbit(irab), rabbit_ages{irab}, sum(num_types)))
	set(gca, 'FontSize', font_size)
	grid on
	box on


end

%%

% Get BMFs/WMFs
BE_MTFs = strcmp(sessions.MTF, 'BE');
BMFs = sessions.BMF(BE_MTFs);
BMFs(isnan(BMFs)) = [];

edges = [0.2 2 4 8 16 32 64 128 254 512 1028];
edges2 = zeros(10, 1);
for iedge = 1:10
	edges2(iedge) = sqrt(edges(iedge)*edges(iedge+1));
end

% BMFs
nexttile()
histogram(BMFs, edges2,'FaceColor', 'k', 'EdgeColor','k')
title('BE: Distribution of BMFs')
hold on
xline(exp(median(log(BMFs(BMFs~=0)))), 'k', 'LineWidth',1.5)
%legend('Distribution', ['Mean=' num2str(round(mean(BMFs))) 'Hz'], ...
%	['Median=' num2str(round(median(BMFs))) 'Hz'])
xlabel('BMF (Hz)')
ylabel('# Neurons')
xticks([2 4 8 16 32 64 128 254 512])
%xticklabels([2 [] 8 [] 32 [] 128 [] 512])
set(gca, 'FontSize', font_size)
set(gca, 'XScale', 'log');

% BMF per rabbits
putative_units = sessions.Putative_Units;
rabbit_units = cellfun(@(p) str2double(p(2:3)), putative_units(1:334), 'UniformOutput', false);
rabbit_putative_units = [rabbit_units{:}];
rabbit = [24, 25, 27, 29];
for irab = 1:4
	rab = rabbit(irab);
	has_rab = rabbit_putative_units==rab & BE_MTFs(1:334)';

	BMFs = sessions.BMF(has_rab);
	BMFs(isnan(BMFs)) = [];

	nexttile
	histogram(BMFs, edges2, 'FaceColor', color{irab}, 'EdgeColor',color{irab})
	hold on
	xline(exp(median(log(BMFs(BMFs~=0)))), 'k', 'LineWidth',1.5)
	xticks([2 4 8 16 32 64 128 254 512])
	%xticklabels([2 [] 8 [] 32 [] 128 [] 512])
	%legend('Distribution', ['Mean=' num2str(round(mean(BMFs))) 'Hz'], ...
	%	['Median=' num2str(round(median(BMFs))) 'Hz'])
	xlabel('BMF (Hz)')
	set(gca, 'FontSize', font_size)
	set(gca, 'XScale', 'log');

end


%%
BS_MTFs = strcmp(sessions.MTF, 'BS');
WMFs = sessions.WMF(BS_MTFs);
WMFs(isnan(WMFs)) = [];


nexttile()
histogram(WMFs, edges2,'FaceColor', 'k', 'EdgeColor','k')
title('BS: Distribution of WMFs')
hold on
xlabel('WMF (Hz)')
xline(exp(median(log(WMFs(WMFs~=0)))), 'k', 'LineWidth',1.5)

%legend('Distribution', ['Mean=' num2str(round(mean(WMFs))) 'Hz'], ...
%	['Median=' num2str(round(median(WMFs))) 'Hz'])
ylabel('# Neurons')
set(gca, 'FontSize', font_size)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');

% WMF per rabbits
putative_units = sessions.Putative_Units;
rabbit_units = cellfun(@(p) str2double(p(2:3)), putative_units(1:334), 'UniformOutput', false);
rabbit_putative_units = [rabbit_units{:}];
rabbit = [24, 25, 27, 29];
for irab = 1:4

	rab = rabbit(irab);
	has_rab = rabbit_putative_units==rab & BS_MTFs(1:334)';
	WMFs = sessions.WMF(has_rab);
	WMFs(isnan(WMFs)) = [];

	nexttile
	histogram(WMFs, edges2, 'FaceColor', color{irab}, 'EdgeColor',color{irab})
	hold on
	xline(exp(median(log(WMFs(WMFs~=0)))), 'k', 'LineWidth',1.5)
	%legend('Distribution', ['Mean=' num2str(round(mean(WMFs))) 'Hz'], ...
	%	['Median=' num2str(round(median(WMFs))) 'Hz'])
	xlabel('WMF (Hz)')
	set(gca, 'FontSize', font_size)
	xticks([2 4 8 16 32 64 128 254 512])
	set(gca, 'XScale', 'log');
end

%% Hybrids

H_MTFs = contains(sessions.MTF, 'H');
WMFs = sessions.WMF(H_MTFs);
WMFs(isnan(WMFs)) = [];
BMFs = sessions.BMF(H_MTFs);
BMFs(isnan(BMFs)) = [];


nexttile()
histogram(BMFs, edges2)
hold on
histogram(WMFs, edges2)
title('Hybrid: BMF & WMF')
xlabel('BMF/WMF (Hz)')
legend('BMF', 'WMF', 'Location','northwest')
ylabel('# Neurons')
set(gca, 'FontSize', font_size)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');

% WMF per rabbits
putative_units = sessions.Putative_Units;
rabbit_units = cellfun(@(p) str2double(p(2:3)), putative_units(1:334), 'UniformOutput', false);
rabbit_putative_units = [rabbit_units{:}];
rabbit = [24, 25, 27, 29];
for irab = 1:4

	rab = rabbit(irab);
	has_rab = rabbit_putative_units==rab & H_MTFs(1:334)';
	WMFs = sessions.WMF(has_rab);
	WMFs(isnan(WMFs)) = [];
	BMFs = sessions.BMF(has_rab);
	BMFs(isnan(BMFs)) = [];

	nexttile
	histogram(BMFs, edges2)
	hold on
	histogram(WMFs, edges2)
	xlabel('BMF/WMF (Hz)')
	set(gca, 'FontSize', font_size)
	xticks([2 4 8 16 32 64 128 254 512])
	set(gca, 'XScale', 'log');
end



%% Get CFs for each putative neuron
clear MTF_type

% WBTIN Diotic
CFs = sessions.CF;

edges = [0 500 1000 2000 4000 8000 13000];
names = categorical({'<500', '500-1k', '1k-2k', '2k-4k', '4k-8k', '8k+'});
names = reordercats(names,{'<500', '500-1k', '1k-2k', '2k-4k', '4k-8k', '8k+'});

CF = CFs;
CF(CF==0) = [];
[N, edges1] = histcounts(CF, edges);

nexttile()
bar(names,N,'FaceColor', 'k', 'EdgeColor','k');
grid on
box on
ylabel('# Neurons')
xlabel('CF (Hz)')
set(gca, 'FontSize', font_size)
title(sprintf('CF Distribution, n=%d', length(CF)))

% CFs per rabbit
for irab = 1:4

	rab = rabbit(irab);
	has_rab = rabbit_putative_units==rab;
	CFs = sessions.CF(has_rab);

	edges = [0 500 1000 2000 4000 8000 13000];
	names = categorical({'<500', '500-1k', '1k-2k', '2k-4k', '4k-8k', '8k+'});
	names = reordercats(names,{'<500', '500-1k', '1k-2k', '2k-4k', '4k-8k', '8k+'});

	CF = CFs;
	CF(CF==0) = [];
	[N, edges1] = histcounts(CF, edges);

	nexttile()
	bar(names,N, 'FaceColor', color{irab}, 'EdgeColor',color{irab});
	grid on
	box on
	xlabel('CF (Hz)')
	set(gca, 'FontSize', font_size)

end

%% Export 

exportgraphics(gcf, fullfile(savepath, 'manuscript', 'data_per_rabbit.png'), 'Resolution', 600)
