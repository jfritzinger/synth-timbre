%% Results3_Comparison
% J. Fritzinger, 1/16/23
clear 

%% Load in spreadsheets 

spreadsheet_name = 'SynthSessions.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

base = getPaths();
fpath = 'data/aro-2024';

%% Find sessions with contra and diotic 

matched_43 = sessions.('43dB_con') & sessions.('43dB');
matched_63 = sessions.('63dB_con') & sessions.('63dB');
matched_83 = sessions.('83dB_con') & sessions.('83dB');

%% Analysis 

for ispl = 1:3
	if ispl == 1
		num_sesh = sum(matched_43);
		index = find(matched_43);

		width_bin43 = NaN(1, num_sesh);
		width_contra43 = NaN(1, num_sesh);
	elseif ispl == 2
		num_sesh = sum(matched_63);
		index = find(matched_63);
		width_bin63 = NaN(1, num_sesh);
		width_contra63 = NaN(1, num_sesh);
	else
		num_sesh = sum(matched_83);
		index = find(matched_83);
		width_bin83 = NaN(1, num_sesh);
		width_contra83 = NaN(1, num_sesh);
	end

	for iclus = 1:num_sesh

		% Load in data 
		session = sessions.Session{index(iclus)};
		TT = sessions.TT(index(iclus));
		N = sessions.N(index(iclus));
		CF = sessions.CF(index(iclus));
		MTF = sessions.MTF(index(iclus));
		filename = sprintf('%s_TT%d_N%d.mat', session, TT, N);
		load(fullfile(base, fpath, filename))
		CFs(iclus) = CF;

		[bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, near_CF]...
		= finddata(params, CF);

		
		MTF_type(iclus) = MTF; 
		if contains(MTF, 'H')
			MTF_type2(iclus) = 3;
		elseif strcmp(MTF, 'BE')
			MTF_type2(iclus) = 1;
		elseif strcmp(MTF, 'BS')
			MTF_type2(iclus) = 2;
		else
			MTF_type2(iclus) = 4;
		end

		if ispl == 1
			ds_bin = find(bin & level_43 & even);
			ds_contra = find(contra & level_43 & even);

			width_bin43(iclus) = data{ds_bin(1)}.width;
			width_contra43(iclus) = data{ds_contra(1)}.width;

			% width_bin43(iclus) = max(data{ds_bin(1)}.rate)/data{ds_bin(1)}.width;
			% width_contra43(iclus) = max(data{ds_contra(1)}.rate)/data{ds_contra(1)}.width;
		elseif ispl == 2
			ds_bin = find(bin & level_63 & even);
			ds_contra = find(contra & level_63 & even);

			width_bin63(iclus) = data{ds_bin(1)}.width;
			width_contra63(iclus) = data{ds_contra(1)}.width;

			% width_bin63(iclus) = max(data{ds_bin(1)}.rate)/data{ds_bin(1)}.width;
			% width_contra63(iclus) = max(data{ds_contra(1)}.rate)/data{ds_contra(1)}.width;
		else
			ds_bin = find(bin & level_83 & even);
			ds_contra = find(contra & level_83 & even);

			width_bin83(iclus) = data{ds_bin(1)}.width;
			width_contra83(iclus) = data{ds_contra(1)}.width;

			if contains(MTF, 'H')
				MTF_type3(iclus) = 3;
			elseif strcmp(MTF, 'BE')
				MTF_type3(iclus) = 1;
			elseif strcmp(MTF, 'BS')
				MTF_type3(iclus) = 2;
			else
				MTF_type3(iclus) = 4;
			end
			% width_bin83(iclus) = max(data{ds_bin(1)}.rate)/data{ds_bin(1)}.width;
			% width_contra83(iclus) = max(data{ds_contra(1)}.rate)/data{ds_contra(1)}.width;
		end


	end
end


%% Scatter plot
num_neurons = length(CFs);

pink = [187, 249, 186]/255;
red = [23, 64, 86]/255;
color_gradient = [linspace(red(1),pink(1),num_neurons)', linspace(red(2),pink(2),num_neurons)',...
	linspace(red(3),pink(3),num_neurons)'];
[~, order] = sort(CFs);

figure('Position',[175,501,875,390])
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact')
for ispl = 1:3

	nexttile
	hold on
	
	% if ispl == 1
	% 	gscatter(width_bin43, width_contra43, MTF_type2, 'filled');
	% 	title('43 dB SPL')
	% elseif ispl == 2
	% 	gscatter(width_bin63, width_contra63, MTF_type2,'filled');
	% 	title('63 dB SPL')
	% else 
	% 	gscatter(width_bin83, width_contra83, MTF_type3,'filled');
	% 	title('83 dB SPL')
	% end

	for ineuron = 1:num_neurons
		if ispl == 1
			scatter(width_bin43(order(ineuron)), width_contra43(order(ineuron)),60,  'filled', 'MarkerFaceColor',color_gradient(ineuron,:));
			title(['43 dB SPL, n=' num2str(length(width_bin43))], 'fontsize', 18)
		elseif ispl == 2
			scatter(width_bin63(order(ineuron)), width_contra63(order(ineuron)),60, 'filled', 'MarkerFaceColor',color_gradient(ineuron,:));
			title(['63 dB SPL, n=' num2str(length(width_bin63))], 'fontsize', 18)
		else
			if order(ineuron) <=51
				scatter(width_bin83(order(ineuron)), width_contra83(order(ineuron)),60, 'filled', 'MarkerFaceColor',color_gradient(ineuron,:));
				title(['83 dB SPL, n=' num2str(length(width_bin83))], 'fontsize', 18)
				%legend
			end
		end
	end

	plot([0 2500], [0 2500], 'k')
	%legend('BE', 'BS', 'H', 'F')
	xlabel('Diotic 1/2 Height BW (Hz)')
	if ispl == 1
		ylabel('Contra 1/2 Height BW (Hz)')
	else
		yticklabels([])
	end

	box on
	xlim([0 2400])
	ylim([0 2400])
	%xlim([0 2])
	%ylim([0 2])
	set(gca, 'fontsize', 16)

end

%% Percentages
spl_names = {'43', '63', '83'};


figure('Position',[175,501,875,200])
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact')
for ispl = 1:3

	nexttile
	hold on
	if ispl == 1
		change = width_contra43-width_bin43;
	elseif ispl == 2
		change = width_contra63-width_bin63;
	else 
		change = width_contra83-width_bin83;
	end

	%edges = linspace(-0.6, 0.6, 17);
	edges = linspace(-1600, 1600, 17);
	histogram(change, edges)
	hold on
	xline(0, 'LineWidth',3)
	xline(mean(change), 'r', 'LineWidth',1.5)
	xline(median(change), '--r', 'LineWidth',1.5)
	 
	if ispl == 1
	legend('', 'No change', 'Mean', 'Median', 'Location','northwest', 'EdgeColor','none', 'fontsize', 12)
	end
	xlabel('Contra - Diotic (Hz)')
	if ispl == 1
		ylabel('# neurons')
	end

	box on
	set(gca, 'fontsize', 16)
	title([spl_names{ispl} ' dB SPL'], 'FontSize',18)

end


%% Significance 
% Perform paired t-test
[h, p, ci, stats] = ttest(width_contra43, width_bin43);

% Display results
fprintf('Paired t-test results:\n');
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
fprintf('Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

% Check for significance
if h
    fprintf('43: The difference is statistically significant.\n\n');
else
    fprintf('43: There is no significant difference.\n\n');
end

% Perform paired t-test
[h, p, ci, stats] = ttest(width_contra63, width_bin63);

% Display results
fprintf('Paired t-test results:\n');
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
fprintf('Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

% Check for significance
if h
    fprintf('63: The difference is statistically significant.\n\n');
else
    fprintf('63: There is no significant difference.\n\n');
end

% Perform paired t-test
[h, p, ci, stats] = ttest(width_contra83, width_bin83);

% Display results
fprintf('Paired t-test results:\n');
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
fprintf('Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

% Check for significance
if h
    fprintf('83: The difference is statistically significant.\n\n');
else
    fprintf('83: There is no significant difference.\n\n');
end

%% Functions

function [bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, near_CF] = finddata(population, CF)

	% Find binaural
	bin = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==2, population);

	% Find contra
	contra = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.binmode==1, population);

	% Find F0 = 100Hz
	F0_100 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==100, population);

	% Find F0 = 200Hz
	F0_200 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && p.Delta_F==200, population);

	% Find 43 db SPL
	level_43 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==43 || p.spl==40), population);

	% Find 63 dB SPL
	level_63 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==63 || p.spl==60), population);

	% Find 73 dB SPL
	level_73 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==73 || p.spl==70), population);

	% Find 83 dB SPL
	level_83 = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && (p.spl==83 || p.spl==80), population);

	% Find even 200Hz fpeak_mid
	even = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);

	% Find nearest to CF 
	datasets = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'SPEC_slide') && ...
		strcmp(p.SPEC_slide_type, 'Spectral_Centroid') && mod(p.fpeak_mid, 200)==0, population);
	near_CF = cellfun(@(p) min(abs(p.fpeak_mid-CF)), population(datasets));
end