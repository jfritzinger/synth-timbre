function distance = dog_objective_function(params, model, Fs, stim, observed_rate, r0, type)

	f = linspace(0, Fs/2, 100000);
    if strcmp(model, 'gaussian')
        fc = 10^params(1);
        sigma = 10^params(2);
        g = params(3);
        W = gaussian_model(f, fc, sigma, g);
    else % DoG model
        W = dog_model(f, params);
    end
    
	nstim = size(stim, 1);
    predicted_rate = zeros(nstim, 1);
    for i = 1:nstim
        predicted_rate(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
    end
    
	% Display parameters 
	% if length(params)==3
	% 	fprintf('Fc = %0.0f, Sigma = %0.2f, g = %0.2f\n', params(1), params(2), params(3))
	% else
	% 	fprintf('g_e=%0.0f, g_i=%0.0f, %c_e=%0.3f, %c_i=%0.3f, f_e=%0.3f, f_i=%0.3f\n',...
		% 	params(1), params(2), 963, params(3), 963, params(4), params(5), params(5))
	% end

	if type == 1 % Distance 
		distance = sum(abs(predicted_rate - observed_rate)); % City-block distance
		%fprintf('Dist = %0.2f\n', distance)
	else % MSE
		%mse = norm(observed_rate-predicted_rate,2)/norm(observed_rate); % Minimize squared error of Euclidian distance (L2)
		distance = 1/length(observed_rate) * sum((observed_rate - predicted_rate).^2); % Compute the mean of the squared errors
		%fprintf('MSE = %0.2f\n', distance)
	end

	% if strcmp(model, 'gaussian')
	% else
	% 	figure
	% 	hold on
	% 	plot(observed_rate)
	% 	plot(predicted_rate)
	% 	x = 1;
	% end

end 