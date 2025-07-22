%% Test Q10

% Human
cf = logspace(log10(700), log10(15000));
q10 = ((cf./1000).^0.3)*12.7*0.505+0.2085;

figure
plot(cf, q10)
set(gca, 'xscale', 'log')
xlim([500 20000])
%ylim([1, 16])
xticklabels([1 10])
xlabel('CF (kHz)')
ylabel('Q10')
hold on

% Human 2
q10 = cf./24.7./(4.37*(cf./1000)+1)*0.505+0.2085;
plot(cf, q10)


% Cat
q10 = 10.^(0.4708*log10(cf./1e3)+0.4664);
plot(cf, q10)

% Rabbit
q10 = (3.62.*(cf./1000).^3.65)./((cf./1000).^3.65+3.27^3.65)+2.022;
%q10 = ((3.62.*(cf./1000).^3.65)./(((cf./1000).^3.65)+75.5265))+2.022;
plot(cf, q10)

% Legend
legend('Human Shera et al. (PNAS 2002)',...
	'Human Glasberg & Moore (Hear. Res. 1990)', 'Cat', 'Rabbit', ...
	'Location','northoutside')
grid on