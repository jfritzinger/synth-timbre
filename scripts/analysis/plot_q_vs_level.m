%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath, "peak_picking.xlsx"));

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

tables.Q(isnan(tables.Q)) = 0;
tables.SPL = categorical(tables.SPL);
tables.SPL = reordercats(tables.SPL, ["43", "63", "73", "83"]);

full_equation = 'Q~SPL+(1+Putative)';
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

%% Plot 

qs_new = qs - qs(:,1);
figure('position', [756,106,586,831])
spls = [43, 63, 73, 83];
differences = NaN(length(qs_new), 1);
for ii = 1:length(qs_new)
	if ~isnan(qs_new(ii,3)) && ~isnan(qs_new(ii,1))
		differences(ii) = qs_new((ii),3)-qs_new((ii),1);
	elseif ~isnan(qs_new(ii,2)) && ~isnan(qs_new(ii,1))
		differences(ii) = qs_new((ii),2)-qs_new((ii),1);
	elseif ~isnan(qs_new(ii,3)) && ~isnan(qs_new(ii,2))
		differences(ii) = qs_new((ii),3)-qs_new((ii),2);
	end
end
decrease = find(sign(differences)==-1);
increase = find(sign(differences)==1);

hold on
yline(0)
plot(spls, qs_new(increase,:)', 'color', [0.85,0.11,0.38, 0.6], 'LineWidth',1.5)
plot(spls, qs_new(decrease,:)', 'color',[0.12,0.53,0.90, 0.6], 'LineWidth',1.5)
xticks(spls)
set(gca, 'fontsize', 18)
xlabel('Level (dB SPL)')
ylabel('Q')
box on
xlim([40 86])
plot(spls, mean(qs_new(increase,:), 'omitnan'), 'k', 'LineWidth',2)
plot(spls, median(qs_new(increase,:), 'omitnan'), '--k', 'LineWidth',2)
plot(spls, mean(qs_new(decrease,:), 'omitnan'), 'k', 'LineWidth',2)
plot(spls, median(qs_new(decrease,:), 'omitnan'), '--k', 'LineWidth',2)

