function W = dog_model(fpeaks, DOGparams)

% Difference of Gaussians (DoG) model

	% excGauss = ge * exp(-(f - fc).^2 / (2 * sigma_e^2));
	% inhGauss = gi * exp(-(f - fc).^2 / (2 * sigma_i^2));
    % W = excGauss - inhGauss;

	% Set Parameters
	s_exc = DOGparams(1);
	s_inh = DOGparams(2);
	sigma_exc = DOGparams(3);
	sigma_inh = DOGparams(4);
	CF_exc = DOGparams(5);
	CF_inh = DOGparams(6);
	gauss_exc = normpdf(fpeaks, CF_exc, sigma_exc);
	gauss_inh = normpdf(fpeaks, CF_inh, sigma_inh);
	gauss_exc = s_exc*(gauss_exc./max(gauss_exc));
	gauss_inh = s_inh*(gauss_inh./max(gauss_inh));
	W = gauss_exc - gauss_inh;

	% Plot to test 
	% figure
	% plot(gauss_exc)
	% hold on
	% plot(-1*gauss_inh)
	% plot(W)
	% xlim([0 10000])
end