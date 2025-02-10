function [gaussian_params, dog_params] = fitGaussAndDoG(params, CF, Fs, observed_rate, r0)
% fitGausAndDoG

type = 2; % 1: distance, 2: MSE
stim = params{1}.stim;
timerVal = tic;

% Fit gaussian
log_CF = log10(CF);
init = [log_CF,  1.4,   100]; % Initial guess (CF, sigma, g)
lb = [log_CF-1,  0.001, 0]; % Lower bounds
ub = [log_CF+1,  4,     Inf]; % Upper bounds

options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
	'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
	'StepTolerance', 1e-16, 'display', 'off');
gaussian_params = fmincon(@(p) ...
	dog_objective_function(p, 'gaussian', Fs, stim, observed_rate, r0, type), ...
    init, [], [], [], [], lb, ub, [], options);


% Fit DoG model
best_fval = Inf;
for istarts = 1:10

	g_exc_init = 100 + (100000 - 100) * rand(1);
	g_inh_init = 100 + (100000 - 100) * rand(1);

	%			g_exc, g_inh, s_exc, s_inh,  CF_exc, CF_inh
	dog_init = [g_exc_init, g_inh_init, 2,     2.5,    log_CF, log_CF]; % Initial guess
	dog_lb = [100,   100,     1,     1,  log_CF-1, log_CF-1]; % Lower bounds
	dog_ub = [100000, 100000, 4,     4,      log_CF+1, log_CF+1]; % Upper bounds

	options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-10, ...
		'MaxFunEvals', 10^10, 'maxiterations', 500, 'ConstraintTolerance', 1e-10, ...
		'StepTolerance', 1e-12, 'display', 'off');
	[dog_params, fval] = fmincon(@(p) dog_objective_function(p, 'dog', Fs, stim, observed_rate, r0, type), ...
		dog_init, [], [], [], [], dog_lb, dog_ub, [], options);
	
	if fval < best_fval
		best_x = dog_params;
		best_fval = fval;
	end
end
dog_params = best_x;
disp(['Model took ' num2str(toc(timerVal)) ' seconds'])




end