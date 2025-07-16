%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
%tables = readtable(fullfile(datapath, "peak_picking2.xlsx"));
tables = readtable(fullfile(datapath,"LMM", "peak_picking_excludeflat.xlsx"));

%% Set up figure

figure('position', [519,299,7*ppi,7*ppi])
tiledlayout(2, 2)
fontsize = 12;
legsize = 10;
labelsize = 24;

%% Scatter plot of CF vs HHBW (all MTF)

spls = [43, 63, 73, 83];
is200 = tables.F0 == 200;
%isMTF = strcmp(tables.MTF, 'BE') | strcmp(tables.MTF, 'BS');
for ibin = 2
	isbin = tables.binmode == ibin;
	for ispl = 2

		% Get data
		islevel = tables.SPL == spls(ispl);
		index = islevel & isbin & is200; % & isMTF;

		% Data
		CFs = tables.CF(index);
		Qs = tables.Q(index);
		MTFs = tables.MTF(index);

		% Plot
		nexttile
		%gscatter(CFs/1000, Qs, MTFs, 'filled')
		scatter(CFs/1000, Qs, 'filled', 'MarkerEdgeColor','k')
		hold on
		box on

		% Fit linear regression line
		mdl = fitlm(log10(CFs), log10(Qs));
		x = 0.3:0.5:10000;
		p(1) = mdl.Coefficients.Estimate(2,1);
		p(2) = mdl.Coefficients.Estimate(1,1);
		p(3) = mdl.Coefficients.pValue(2);
		p(4) = mdl.Rsquared.Ordinary;
		mdlfit(ibin, ispl,:) = 10.^(p(1)*log10(x)+p(2));
		mdlplot = squeeze(mdlfit(ibin, ispl, :));
		plot(x/1000, mdlplot, 'k');

		% Plot labels 
		number = Qs;
		number(isnan(number)) = [];
		xlabel('CF (kHz)')
		ylabel('Q')
		ylim([0.35 50])
		xlim([0.3 10])
		set(gca, 'XScale', 'log', 'YScale', 'log')
		xticks([0 200 500 1000 2000 5000 10000]/1000)
		yticks([0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000])
		grid on
		set(gca, 'fontsize', fontsize)
		msg = ['n=' num2str(length(number))];
		text(0.05, 0.95, msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
	end
end
box off

%% Significant interactions between level and CF grouping 

is200 = tables.F0 == 200;
spls = [43, 63, 73, 83];
isbin = tables.binmode == 2; % | tables.binmode == 1;
CFgroup = {'Low', 'Med', 'High'};
for iCF = 1:3
	for ispl = 1:4

		% Get data
		islevel = tables.SPL == spls(ispl);
		isCFgroup = strcmp(CFgroup{iCF}, tables.CF_Group);
		index = islevel & isbin & is200 & isCFgroup;

		% Data
		CFs = tables.CF(index);
		Qs = tables.Q(index);
		Q(iCF, ispl) = mean(Qs);
		Q_sem(iCF, ispl) = std(Qs)/sqrt(length(Qs));
	end
end

nexttile
errorbar(Q', Q_sem', 'LineWidth',2)
%xlabel('CF Group')
xlabel('Level (dB SPL)')
xlim([0.5 4.5])
xticks(1:4)
%xticklabels({'Low', 'Med', 'High'})
xticklabels([43, 63, 73, 83])
ylabel('Q')
ylim([0 11])
%legend('43 dB SPL', '63 dB SPL', '73 dB SPL', '83 dB SPL', 'Location','best')
legend('Low CF', 'Med CF', 'High CF', 'Location','best', 'fontsize',legsize)
set(gca, 'fontsize', fontsize)
grid on
box off

%% Plot q vs MTF 

ispeak = strcmp(tables.Type, 'Peak');
isdip = strcmp(tables.Type, 'Dip');
isflat = strcmp(tables.Type, 'Flat');
signed_Q = tables.Q;
signed_Q(isdip) = signed_Q(isdip) * -1;
signed_Q(isflat) = 0;
signed_MTF = tables.MTF_str;

isbin = tables.binmode == 2;
is200 = tables.F0 == 200;
spls = [43, 63, 73, 83];
for ispl = 2

	% Get data
	islevel = tables.SPL == spls(ispl);
	index = islevel & isbin & is200;
	isBE = strcmp(tables.MTF(index), 'BE');
	isBS = strcmp(tables.MTF(index), 'BS');
	isH = contains(tables.MTF(index), 'H');
	
	% Data
	CFs = tables.CF(index);
	Qs = tables.Q(index);
	MTFs = tables.MTF(index);
	MTF_str = signed_MTF(index);

	% Plot
	nexttile
	hold on
	scatter(MTF_str(isBS), Qs(isBS), 'filled', 'MarkerEdgeColor','k')
	scatter(MTF_str(isBE), Qs(isBE), 'filled', 'MarkerEdgeColor','k')
	%scatter(MTF_str(isH), Qs(isH), 'filled', 'MarkerEdgeColor','k')
	xline(0)
	yline(0)
	grid on
	ylabel('Q')
	xlabel('MTF % Change')
	legend(['BS, n=' num2str(length(Qs(isBS)))],...
		['BE, n=' num2str(length(Qs(isBE)))],'fontsize', legsize)
	set(gca, 'fontsize', fontsize)
end
box off


%% Significant interaction between MTF and binmode 

is200 = tables.F0 == 200;
isBE = strcmp(tables.MTF, 'BE');
isBS = strcmp(tables.MTF, 'BS');
isH = contains(tables.MTF, 'H');
isF = strcmp(tables.MTF, 'F');
for ibin = 1:2

	% Get data
	isbin = tables.binmode == ibin;
	ind_BE = isbin & is200 & isBE;
	ind_BS = isbin & is200 & isBS;
	ind_H = isbin & is200 & isH;
	ind_F = isbin & is200 & isF;

	% Data
	Q_all2(ibin,:) = [mean(tables.Q(ind_BE)) mean(tables.Q(ind_BS)) ...
		mean(tables.Q(ind_H)) mean(tables.Q(ind_F))];
	Q_sem2(ibin,:) = [std(tables.Q(ind_BE))/sqrt(length(tables.Q(ind_BE)))...
		std(tables.Q(ind_BS))/sqrt(length(tables.Q(ind_BS)))...
		std(tables.Q(ind_H))/sqrt(length(tables.Q(ind_H)))...
		std(tables.Q(ind_F))/sqrt(length(tables.Q(ind_F)))];
end

% Plot
nexttile
hold on
errorbar(Q_all2', Q_sem2', 'LineWidth',2)
xticks(1:4)
xticklabels({'BE', 'BS', 'Hybrid', 'Flat'})
xlim([0.5 4.5])
ylim([0 6.5])
ylabel('Q')
xlabel('MTF Type')
legend('Contra', 'Diotic', 'Location','best', 'fontsize', legsize)
set(gca, 'fontsize', fontsize)
grid on

%% Significant interaction between MTF and binmode 
% 
% figure('position', [519,299,1150,250])
% tiledlayout(1, 4, 'TileSpacing','tight', 'Padding','compact')
% 
% isBE = strcmp(tables.MTF, 'BE');
% isBS = strcmp(tables.MTF, 'BS');
% is200 = tables.F0 == 200;
% spls = [43, 63, 73, 83];
% isbin = tables.binmode == 2;
% iscontra = tables.binmode == 1;
% for ispl = 1:4
% 
% 	% Get data
% 	islevel = tables.SPL == spls(ispl);
% 	index = islevel & is200;
% 
% 	% Data
% 	Qs_BE_contra = tables.Q(index & isBE & iscontra);
% 	Qs_BE_bin = tables.Q(index & isBE & isbin);
% 	putative_BE_contra = tables.Putative(index & isBE & iscontra);
% 	putative_BE_bin = tables.Putative(index & isBE & isbin);
% 	[v_bin, v_contra] = matchBinAndContra(putative_BE_bin, putative_BE_contra);
% 
% 	Qs_BS_contra = tables.Q(index & isBS & iscontra);
% 	Qs_BS_bin = tables.Q(index & isBS & isbin);
% 	putative_BS_contra = tables.Putative(index & isBS & iscontra);
% 	putative_BS_bin = tables.Putative(index & isBS & isbin);
% 	[v_bin_BS, v_contra_BS] = matchBinAndContra(putative_BS_bin, putative_BS_contra);
% 
% 	% Plot
% 	nexttile
% 	hold on
% 	scatter(Qs_BE_contra(v_contra), Qs_BE_bin(v_bin), 'filled', 'MarkerEdgeColor','k')
% 	scatter(Qs_BS_contra(v_contra_BS), Qs_BS_bin(v_bin_BS), 'filled', 'MarkerEdgeColor','k')
% 	xline(0)
% 	yline(0)
% 	box on
% 	ylabel('Q')
% 	xlabel('MTF % Change')
% 	legend('BS', 'BE')
% 	set(gca, 'fontsize', fontsize)
% 
% end

%% Annotations and labels 

% Set annotations
left = [0.03 0.51 0.03 0.51];
bottom = [0.95 0.95 0.47 0.47];
label = {'A', 'B', 'C', 'D'};
for ii = 1:4
	annotation('textbox',[left(ii) bottom(ii) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end

%% FUNCTIONS 

function [v_bin, v_contra] = matchBinAndContra(putative_bin, putative_contra)

v_bin = zeros(1, length(putative_bin));
v_contra = zeros(1, length(putative_contra));
for ii = 1:length(putative_contra)
	putative = putative_contra{ii};
	for jj = 1:length(putative_bin)
		matched = strcmp(putative, putative_bin{jj});
		if matched == 1
			v_contra(ii) = 1;
			v_bin(jj) = 1;
		end
	end
end
v_bin = logical(v_bin);
v_contra = logical(v_contra);

end
