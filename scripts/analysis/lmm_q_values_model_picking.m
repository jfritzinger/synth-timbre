%% lmm_q_values_model_picking.m
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
tables_old = readtable(fullfile(datapath, "peak_picking.xlsx"));
%tables_old.Q_log = log10(tables_old.Q);

% Deal with NaNs
%tables.Q(isnan(tables.Q)) = 0;
index = isnan(tables_old.Q);
tables = tables_old(~index, :);

tables.SPL = categorical(tables.SPL);
tables.SPL = reordercats(tables.SPL, ["43", "63", "73", "83"]);



%% Check normality

% Not normal 
figure('Position',[560,526,894,322])
tiledlayout(1, 2)

nexttile
histogram(tables.Q)
title('Histogram of Q values')
xlabel('Q values')
ylabel('# Units')
set(gca, 'fontsize', 14)

% Normalize 
nexttile
histogram(tables.Q_log)
%tables.Q_log = log10(tables.Q);
title('Histogram of log-transformed Q values')
xlabel('log10(Q values)')
ylabel('# Units')
set(gca, 'fontsize', 14)

%% Set up equations

% Only random effect
equation{1} = '~1+CF_Group+(1|Putative)+(CF_Group|Putative)';

% One fixed effect 
equation{2} = '~SPL+(1|Putative)';
equation{3} = '~MTF+(1|Putative)';
equation{4} = '~F0+(1|Putative)';
equation{5} = '~CF_Group+(1|Putative)';
equation{6} = '~binmode+(1|Putative)';

% Two fixed effects 
equation{7} = '~SPL+MTF+(1|Putative)';
equation{8} = '~SPL+CF_Group+(1|Putative)';

% All fixed effects 
equation{9} = '~SPL*CF_Group+MTF+binmode+F0+(1|Putative)';
equation{10} = '~SPL*CF_Group*binmode+MTF+F0+(1|Putative)';
equation{11} = '~SPL*CF_Group*binmode+MTF+F0+(1|Putative)';
equation{12} = '~SPL*CF_Group*binmode+MTF*SPL+(1|Putative)';
equation{13} = '~SPL*CF_Group*binmode+MTF*SPL+F0*CF_Group+(1|Putative)';
equation{14} = '~SPL*CF_Group*binmode+MTF*SPL+MTF*CF_Group+F0*CF_Group+(1|Putative)';


% JUST stimulus parameters 
equation{15} = '~binmode*F0*SPL+(1|Putative)';


%% Test 

ref =  1; % 9, 13, 16
test = 5; 

full_equation = ['Q_log' equation{ref}];
mdl = fitlme(tables, full_equation,'FitMethod','ML', 'DummyVarCoding','effects');
R2 = mdl.Rsquared.Ordinary;

full_equation = ['Q_log' equation{test}];
altmdl = fitlme(tables, full_equation,'FitMethod','ML', 'DummyVarCoding','effects');
altR2 = altmdl.Rsquared.Ordinary;

disp('----------------------------------------------------------------------------------------------------------')
a = compare(mdl, altmdl);
disp(a)
fprintf('Model 1: R^2 = %0.6f\nModel 2: R^2 = %0.6f\n', R2, altR2)
disp('----------------------------------------------------------------------------------------------------------')

%% 

full_equation = ['Q' equation{2}]; 
mdl = fitlme(tables, full_equation,'FitMethod','REML', 'DummyVarCoding','effects');
R2 = mdl.Rsquared.Ordinary;
disp(mdl)
anova(mdl, 'DFMethod','satterthwaite')


figure('Position',[276,526,1048,322])
tiledlayout(1, 3)

nexttile
plotResiduals(mdl, 'fitted')

nexttile
plotResiduals(mdl,'probability')

F = fitted(mdl);
R = response(mdl);
nexttile 
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')