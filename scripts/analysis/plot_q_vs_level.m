%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath,"LMM", "peak_picking_excludeflat.xlsx"));

%% Change in Q vs level 

% Find sessions for target synthetic timbre response
all_neurons = tables.Putative;
neurons = unique(all_neurons);
num_units = size(neurons, 1);
isbin = tables.binmode == 2;
is200 = tables.F0 == 200;

SPLs = [43, 63, 73, 83];
qs = NaN(num_units, 4);
for isesh = 1:num_units

	putative = neurons{isesh};
	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);

	for ispl = 1:4
		ind = isput & isbin & is200 & tables.SPL==SPLs(ispl);
		if any(ind)
			qs(isesh, ispl) = tables.Q(ind);
			qs_log(isesh, ispl) = tables.Q_log(ind);
			CF_group(isesh) = tables.CF_Group(ind);
		end
	end
end


%% Plot box chart & anova 

[~, ~, stats] = anova1(qs);
multcompare(stats)

figure
boxplot(qs)

figure
plot(mean(qs, 'omitnan'))

figure
anova1(qs)

%% 

% tables.SPL = categorical(tables.SPL);
% tables.SPL = reordercats(tables.SPL, ["43", "63", "73", "83"]);

full_equation = 'Q_log~SPL+(1+Putative)';
mdl = fitlme(tables, full_equation,'FitMethod','REML', 'DummyVarCoding','effects');
R2 = mdl.Rsquared.Ordinary;
disp(mdl)
anova(mdl, 'DFMethod','satterthwaite')

figure
plotResiduals(mdl, 'fitted')
figure
plotResiduals(mdl,'probability')

F = fitted(mdl);
R = response(mdl);
figure();
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')

%% Histograms 

edges = linspace(0, 15, 20);
figure
histogram(qs(:,1), edges)
hold on
histogram(qs(:,4), edges)

figure
histogram(qs(:,1)-qs(:,4))

%% Older plots 

figure('position', [890,461,775,305])
tiledlayout(1, 2, 'Padding','compact', 'TileSpacing','compact')
spls = [43, 63, 73, 83];
differences = NaN(length(qs), 1);
for ii = 1:length(qs)
	if ~isnan(qs(ii,3)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),3)-qs((ii),1);
	elseif ~isnan(qs(ii,2)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),2)-qs((ii),1);
	elseif ~isnan(qs(ii,3)) && ~isnan(qs(ii,2))
		differences(ii) = qs((ii),3)-qs((ii),2);
	end
end
decrease = find(sign(differences)==-1);
increase = find(sign(differences)==1);

nexttile
hold on
plot(spls, qs(increase,:)', 'color', [0.85,0.11,0.38, 0.6], 'LineWidth',1.5)
xticks(spls)
set(gca, 'fontsize', 18)
title(['Sharpening Peak (n=' num2str(length(increase)) ')'])
xlabel('Level (dB SPL)')
ylabel('Q')
box on
xlim([40 86])
plot(spls, mean(qs(increase,:), 'omitnan'), 'k', 'LineWidth',2)
plot(spls, median(qs(increase,:), 'omitnan'), '--k', 'LineWidth',2)
% [~, ~, stats] = anova1(qs(increase,:));
% multcompare(stats)

nexttile
hold on
plot(spls, qs(decrease,:)', 'color',[0.12,0.53,0.90, 0.6], 'LineWidth',1.5)
xticks(spls)
set(gca, 'fontsize', 18)
title(['Broadening Peak (n=' num2str(length(decrease)) ')'])
xlabel('Level (dB SPL)')
box on
xlim([40 86])
yticklabels([])
plot(spls, mean(qs(decrease,:), 'omitnan'), 'k', 'LineWidth',2)
plot(spls, median(qs(decrease,:), 'omitnan'), '--k', 'LineWidth',2)
legend('', '', '', '', '', '', '', '', '', '','','','','','','', '', '','Mean', 'Median', 'fontsize', 18)
% [~, ~, stats] = anova1(qs(decrease,:));
% multcompare(stats)

%% CURRENT 

ilow = cellfun(@(d) strcmp('Low', d), CF_group);
imed = cellfun(@(d) strcmp('Med', d), CF_group);
ihigh = cellfun(@(d) strcmp('High', d), CF_group);

spls = [43, 63, 73, 83];
differences = NaN(length(qs), 1);
for ii = 1:length(qs)
	if ~isnan(qs(ii,3)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),3)-qs((ii),1);
	elseif ~isnan(qs(ii,2)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),2)-qs((ii),1);
	elseif ~isnan(qs(ii,3)) && ~isnan(qs(ii,2))
		differences(ii) = qs((ii),3)-qs((ii),2);
	end
end
decrease = sign(differences)==-1;
increase = sign(differences)==1;

figure
tiledlayout(3, 2, 'TileSpacing','compact')
for igroup = 1:3
	
	if igroup == 1
		inc = qs(ilow' & increase,:);
		dec = qs(ilow' & decrease,:);
	elseif igroup == 2
		inc = qs(imed' & increase,:);
		dec = qs(imed' & decrease,:);
	else
		inc = qs(ihigh' & increase,:);
		dec = qs(ihigh' & decrease,:);
	end

	nexttile
	hold on
	for i = 1:4
		swarmchart(ones(length(inc),i)*i, inc(:,i), 15, 'filled')
	end
	boxplot(inc)
	ylim([0 15])
	ylabel('Q-value')
	xticklabels([43, 63, 73, 83])
	xlabel('Sound Level (dB SPL)')
	title('Sharpening Peak')
	set(gca, 'fontsize', 16)
	plot(1:4, median(inc, 'omitnan'), '--k', 'LineWidth',1.5)

	nexttile
	hold on
	for i = 1:4
		swarmchart(ones(length(dec),i)*i, dec(:,i), 15, 'filled')
	end
	boxplot(dec)
	ylim([0 15])
	xticklabels([43, 63, 73, 83])
	xlabel('Sound Level (dB SPL)')
	title('Broadening Peak')
	set(gca, 'fontsize', 16)
	plot(1:4, median(dec, 'omitnan'), '--k', 'LineWidth',1.5)
end

%% Plot 
% 
% qs_new_log = qs_log - qs_log(:,1);
% qs_new = qs - qs(:,1);
% figure('position', [756,106,586,831])
% spls = [43, 63, 73, 83];
% differences2 = NaN(length(qs_new), 1);
% for ii = 1:length(qs_new)
% 	if ~isnan(qs_new(ii,1))
% 		diff_temp = qs_new(ii,:) - qs_new(ii,1);
% 	elseif ~isnan(qs_new(ii,2))
% 		diff_temp = qs_new(ii,:) - qs_new(ii,2);
% 	elseif ~isnan(qs_new(ii,3))
% 		diff_temp = qs_new(ii,:) - qs_new(ii,3);
% 	end
% 	[max_diff, max_index] = max(abs(diff_temp));
% 	differences2(ii) = diff_temp(max_index);
% end
% decrease = find(sign(differences2)==-1);
% increase = find(sign(differences2)==1);
% 
% hold on
% yline(0)
% plot(spls, qs_new(increase,:)', 'color', [0.85,0.11,0.38, 0.6], 'LineWidth',1.5)
% plot(spls, qs_new(decrease,:)', 'color',[0.12,0.53,0.90, 0.6], 'LineWidth',1.5)
% xticks(spls)
% set(gca, 'fontsize', 18)
% xlabel('Level (dB SPL)')
% ylabel('Q')
% box on
% xlim([40 86])
% plot(spls, mean(qs_new(increase,:), 'omitnan'), 'k', 'LineWidth',2)
% plot(spls, median(qs_new(increase,:), 'omitnan'), '--k', 'LineWidth',2)
% plot(spls, mean(qs_new(decrease,:), 'omitnan'), 'k', 'LineWidth',2)
% plot(spls, median(qs_new(decrease,:), 'omitnan'), '--k', 'LineWidth',2)

%% Linear regression for each neuron to check significance 
% 
% x = 1:4;
% num_neurons = size(qs, 1);
% p_values = zeros(num_neurons, 1);
% slopes = zeros(num_neurons, 1);
% for i = 1:num_neurons
% 
% 	y = qs_log(i, :)';
%     tbl = table(x', y, 'VariableNames', {'X', 'Q'});
%     %mdl = fitlm(tbl, 'Q ~ X');
% 	mdl = fitlm(tbl,'Q ~ X + X^2');
% 
%     p_values(i) = mdl.Coefficients.pValue(2); % p-value
%     slopes(i) = mdl.Coefficients.Estimate(2); % slope
% 
% 	% Plot
% 	% figure
% 	% plot(x, y)
% 	% hold on
% 	% y_fit = predict(mdl, table(x', 'VariableNames', {'X'}));
% 	% plot(x, y_fit, 'r-', 'LineWidth', 2);
% end
% 
% alpha = 0.05; % Significance level
% significant_neurons = sum(p_values < alpha);
% 
% increasing_neurons = sum(slopes > 0 & p_values < alpha);
% decreasing_neurons = sum(slopes < 0 & p_values < alpha);
% 
% fprintf('Total neurons with significant changes: %d\n', significant_neurons);
% fprintf('Neurons with significant increase: %d\n', increasing_neurons);
% fprintf('Neurons with significant decrease: %d\n', decreasing_neurons);

%% Another plot 
% 
% 
% figure('position', [890,461,775,305])
% tiledlayout(1, 2, 'Padding','compact', 'TileSpacing','compact')
% spls = [43, 63, 73, 83];
% 
% nexttile
% hold on
% plot(spls, qs_log', 'color', [0.85,0.11,0.38, 0.6], 'LineWidth',1.5)
% xticks(spls)
% set(gca, 'fontsize', 18)
% xlabel('Level (dB SPL)')
% ylabel('Q')
% box on
% xlim([40 86])
% plot(spls, mean(qs(increase,:), 'omitnan'), 'k', 'LineWidth',2)
% plot(spls, median(qs(increase,:), 'omitnan'), '--k', 'LineWidth',2)
% [~, ~, stats] = anova1(qs(increase,:));
% multcompare(stats)

