function [gaussian_params, dog_params] = fitGaussAndDoG(params, CF, Fs, observed_rate, r0)
% fitGausAndDoG

type = 2; % 1: distance, 2: MSE
stim = params{1}.stim;
timerVal = tic;

% Fit gaussian
best_fval = Inf;
for istarts = 1
	log_CF = log10(CF);
	s_init = 1 + (4 - 1) * rand(1);
	g_init = 1000*rand(1);
	init = [log_CF,  s_init,  g_init]; % Initial guess (CF, sigma, g)
	lb = [log_CF-1,  1, 0]; % Lower bounds
	ub = [log_CF+1,  4,     Inf]; % Upper bounds

	options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-10, ...
		'MaxFunEvals', 10^10, 'maxiterations', 500, 'ConstraintTolerance', 1e-10, ...
		'StepTolerance', 1e-10, 'display', 'off');
	[gaussian_params, fval] = fmincon(@(p) ...
		dog_objective_function(p, 'gaussian', Fs, stim, observed_rate, r0, type), ...
		init, [], [], [], [], lb, ub, [], options);

	if fval < best_fval
		best_x = gaussian_params;
		best_fval = fval;
	end
end
gaussian_params = best_x;

% Fit DoG model
best_fval = Inf;
for istarts = 1:50

	g_exc_init = 100 + (100000 - 100) * rand(1);
	g_inh_init = 100 + (100000 - 100) * rand(1);
	s_exc_nit = 1 + (4 - 1) * rand(1);
	s_inh_init = 1 + (4 - 1) * rand(1);

	%			g_exc, g_inh, s_exc, s_inh,  CF_exc, CF_inh
	dog_init = [g_exc_init, g_inh_init, s_exc_nit, s_inh_init, log_CF, log_CF]; % Initial guess
	dog_lb = [100,   100,     1,     1,  log_CF-1, log_CF-1]; % Lower bounds
	dog_ub = [100000, 100000, 4,     4,      log_CF+1, log_CF+1]; % Upper bounds

	options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-15, ...
		'MaxFunEvals', 10^15, 'maxiterations', 800, 'ConstraintTolerance', 1e-15, ...
		'StepTolerance', 1e-15, 'display', 'off');
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