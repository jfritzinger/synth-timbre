%% quantify_peaks.m
%
% Script that plots pie charts and histograms of the peak/dip/flat
% quantifications.
%
%
% Author: J. Fritzinger
% Created: ------; Last revision: 2024-10-22
%
% -------------------------------------------------------------------------
clear

%% Load in data

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, 'peak_picking.xlsx'));

%% Plot all as histogram 

% Plot histogram responses
figure('Position',[68,790,890,220])
tiledlayout(1, 4)
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = tables.binmode == 2;
is200 = tables.F0 == 200;
clear percent_peak percent_dip percent_flat
for ispl = 1:4
	isSPL = tables.SPL == spl(ispl);

	index = isSPL & isBin & is200;

	num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
	num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
	num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));
	all = sum([num_peak num_dip num_flat]);

	percent_peak = num_peak;%/all*100;
	percent_dip = num_dip;%/all*100;
	percent_flat = num_flat;%/all*100;
	
	percent_all = [percent_peak; percent_dip; percent_flat]';

	nexttile
	bar(percent_all, 'stacked')
	title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels({'Peak', 'Dip', 'Flat'})
	ylabel('# Neurons')
	xlabel('Level (dB SPL)')
	%ylim([0 80])

end

%% 

% Plot histogram responses
figure()
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = tables.binmode == 2;
clear percent_peak percent_dip percent_flat
for ispl = 1:4
	isSPL = tables.SPL == spl(ispl);

	index = isSPL & isBin;

	num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
	num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
	num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));
	all = sum([num_peak num_dip num_flat]);

	percent_peak(ispl) = num_peak;%/all*100;
	percent_dip(ispl) = num_dip;%/all*100;
	percent_flat(ispl) = num_flat;%/all*100;


end
percent_all = [percent_peak; percent_dip; percent_flat]';

bar(percent_all, 'stacked')
title([num2str(spl(ispl)) ' dB SPL'])
xticklabels({'43', '63', '73', '83'})
legend({'Peak', 'Dip', 'Flat'}, 'Location','southeast')
ylabel('# Neurons')
xlabel('Level (dB SPL)')
%ylim([0 80])

%% Pie charts

figure('Position',[509,581,731,297])
tiledlayout(2, 4)
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
for iMTF = 1:2

	if iMTF == 1
		MTF_target = 'BS';
	else
		MTF_target = 'BE';
	end
	isMTF = strcmp(tables.MTF, MTF_target);

	for ispl = 1:4
		isSPL = tables.SPL == spl(ispl);
		index = isSPL & isMTF;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));

		nexttile
		piechart([num_peak num_dip num_flat])
		title([MTF_target num2str(spl(ispl)) ' dB SPL'])
		%piechart(tables(index),'Type')
		%legend('Peak', 'Dip', 'Flat')
		%legend

	end
end


%% Histograms

figure('Position',[68,790,890,220])
tiledlayout(1, 4)
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = tables.binmode == 2;
for ispl = 2
	isSPL = tables.SPL == spl(ispl);

	for iMTF = 1:4

		if iMTF == 1
			MTF_target = 'BS';
			isMTF = strcmp(tables.MTF, MTF_target);
		elseif iMTF == 2
			MTF_target = 'BE';
			isMTF = strcmp(tables.MTF, MTF_target);
		elseif iMTF == 3
			MTF_target = 'Hybrid';
			isMTF = contains(tables.MTF, 'H');
		else
			MTF_target = 'F';
			isMTF = strcmp(tables.MTF, MTF_target);
		end
		index = isSPL & isMTF & isBin;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));
		all = sum([num_peak num_dip num_flat]);

		percent_peak(iMTF) = num_peak/all*100;
		percent_dip(iMTF) = num_dip/all*100;
		percent_flat(iMTF) = num_flat/all*100;
	end
	percent_all = [percent_peak; percent_dip; percent_flat]';

	nexttile
	bar(percent_all, 'stacked')
	title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels({'BS', 'BE', 'Hybrid', 'Flat'})
	legend('Peak', 'Dip', 'Flat', 'Location','south')
	ylabel('% Neurons')
	xlabel('MTF Type')
	ylim([0 100])

end

%% # of neurons for BE, BS, flat at 200Hz modulation
% 
% for iMTF = 1:3
% 
% 	if iMTF == 1
% 		MTF_target = 'BS';
% 	elseif iMTF == 2
% 		MTF_target = 'BE';
% 	else
% 		MTF_target = 'F';
% 	end
% 
% 	isMTF = strcmp(tables.MTF_at200, MTF_target);
% 	disp([MTF_target ': ' num2str(sum(isMTF))])
% end

%% Histograms for responses at 200Hz modulation

figure('Position',[509,581,731,297])
tiledlayout(1, 4)
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = tables.binmode == 2;
for iMTF = 1:4

	if iMTF == 1
		MTF_target = 'BS';
	elseif iMTF == 2
		MTF_target = 'BE';
	else
		MTF_target = 'F';
	end
	isMTF = strcmp(tables.MTF_at200, MTF_target);

	for ispl = 1:4
		isSPL = tables.SPL == spl(ispl);
		index = isSPL & isMTF & isBin;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));
		all = sum([num_peak num_dip num_flat]);

		percent_peak(ispl) = num_peak/all*100;
		percent_dip(ispl) = num_dip/all*100;
		percent_flat(ispl) = num_flat/all*100;
	end
	percent_all = [percent_peak; percent_dip; percent_flat]';

	nexttile
	bar(percent_all, 'stacked')
	title(MTF_target)
	xticklabels([43, 63, 73, 83])
	legend('Peak', 'Dip', 'Flat')
	ylabel('Percentage of Neurons')
	xlabel('Level (dB SPL)')

end


%% Histograms for CF ranges: low, medium, and high

% Loop through to find low, medium, and high CFs
CFs = tables.CF;
ilow = CFs<2000;
imed = CFs>=2000 & CFs<4000;
ihigh = CFs>=4000;

% Plot histogram responses
figure('Position',[68,790,890,220])
tiledlayout(1, 4)
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
clear percent_peak percent_dip percent_flat
for ispl = 1:4
	isSPL = tables.SPL == spl(ispl);
	for iCF = 1:3

		if iCF == 1
			CF_target = 'Low';
			isCF = ilow;
		elseif iCF == 2
			CF_target = 'Medium';
			isCF = imed;
		else
			CF_target = 'High';
			isCF = ihigh;
		end

		index = isSPL & isCF & isBin;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), tables.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), tables.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), tables.Type(index)));
		all = sum([num_peak num_dip num_flat]);

		percent_peak(iCF) = num_peak%/all*100;
		percent_dip(iCF) = num_dip%/all*100;
		percent_flat(iCF) = num_flat%/all*100;
	end
	percent_all = [percent_peak; percent_dip; percent_flat]';

	nexttile
	bar(percent_all, 'stacked')
	title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels({'Low', 'Med', 'High'})
	legend({'Peak', 'Dip', 'Flat'}, 'Location','southeast')
	ylabel('# Neurons')
	xlabel('Level (dB SPL)')
	ylim([0 80])

end

%% Violin plots (except need to update matlab first lol)


figure('Position',[328,101,1084,846])
tiledlayout(4, 4, "TileSpacing","tight")
spl = [43, 63, 73, 83];

for ispl = 1:4
	isSPL = tables.SPL == spl(ispl);
	index = isSPL;

	
	MTF_types = tables.MTF(index);
	CFs = tables.CF(index);
	types = tables.Type(index);

	for iMTF = 1:4
		nexttile
		if iMTF == 1
			ind_MTF = strcmp(MTF_types, 'BS');
		elseif iMTF == 2
			ind_MTF = strcmp(MTF_types, 'BE');
		elseif iMTF == 3
			ind_MTF = contains(MTF_types, 'H');
		else
			ind_MTF = strcmp(MTF_types, 'F');
		end

		for itype = 1:3

			if itype == 1
				ind = strcmp(types, 'Peak');
			elseif itype == 2
				ind = strcmp(types, 'Dip');
			else
				ind = strcmp(types, 'Flat');
			end
			CFs_new = CFs(ind&ind_MTF);
			x = itype * ones(1, length(CFs_new));

			hold on
			swarmchart(x, CFs_new./1000, 'filled', 'MarkerEdgeColor','k')
			ylabel('CF (kHz)')
			xlabel('Quantification')
			xticks([1, 2, 3])
			xticklabels({'Peak', 'Dip', 'Flat'})
			yline(2, ':')
			yline(4, ':')
			set(gca, 'yscale', 'log')
			ylim([0.3 10])
			yticks([0.1 0.2 0.5 1 2 5 10])
			
		end
	end
end

% Titles
titles_x = {'BS', 'BE', 'Hybrid', 'Flat'};
locs = linspace(0.16, 0.772, 4);
for ii = 1:4
	annotation('textbox',[locs(ii) 0.933 0.126 0.0459],...
		'String',titles_x{ii},'FontWeight','bold',...
		'FontSize',20,'EdgeColor','none');
end

titles_y = {'43', '63', '73', '83'};
locs = linspace(0.17, 0.82, 4);
for ii = 1:4
	annotation('textbox',[0.02 locs(ii) 0.126 0.0459],...
		'String',titles_y{ii},'FontWeight','bold',...
		'FontSize',20,'EdgeColor','none');
end


%% 
% 
% figure('Position',[328,101,1084,846])
% tiledlayout(4, 3, "TileSpacing","tight")
% spl = [43, 63, 73, 83];
% 
% for ispl = 1:4
% 	isSPL = tables.SPL == spl(ispl);
% 	index = isSPL;
% 
% 
% 	MTF_types = tables.MTF(index);
% 	CFs = tables.CF(index);
% 	types = tables.Type(index);
% 
% 	for itype = 1:3
% 
% 		nexttile
% 		if itype == 1
% 			ind = strcmp(types, 'Peak');
% 		elseif itype == 2
% 			ind = strcmp(types, 'Dip');
% 		else
% 			ind = strcmp(types, 'Flat');
% 		end
% 
% 		for iMTF = 1:4
% 
% 			if iMTF == 1
% 				ind_MTF = strcmp(MTF_types, 'BS');
% 			elseif iMTF == 2
% 				ind_MTF = strcmp(MTF_types, 'BE');
% 			elseif iMTF == 3
% 				ind_MTF = contains(MTF_types, 'H');
% 			else
% 				ind_MTF = strcmp(MTF_types, 'F');
% 			end
% 			CFs_new = CFs(ind&ind_MTF);
% 			x = iMTF * ones(1, length(CFs_new));
% 
% 			hold on
% 			swarmchart(x, CFs_new./1000, 'filled', 'MarkerEdgeColor','k')
% 			ylabel('CF (kHz)')
% 			xlabel('Quantification')
% 			xticks([1, 2, 3, 4])
% 			xticklabels({'BS', 'BE', 'Hybrid', 'Flat'})
% 			yline(2, ':')
% 			yline(4, ':')
% 			set(gca, 'yscale', 'log')
% 			ylim([0.3 10])
% 			yticks([0.1 0.2 0.5 1 2 5 10])
% 			if itype == 1
% 				title('Peak')
% 			elseif itype == 2
% 				title('Dip')
% 			elseif itype == 3
% 				title('Flat')
% 			end
% 
% 		end
% 	end
% end
