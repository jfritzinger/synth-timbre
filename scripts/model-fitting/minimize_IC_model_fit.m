function mse = minimize_IC_model_fit(data, AN, stim_params, model_params, x)

% Run model 
CS_params = [x(1:2) 0.001];
BMFs = x(3:5);
nstim = size(stim_params, 2);
model_outputs = cell(nstim, 1);

parfor istim = 1:nstim
	param = stim_params{istim};
	an_sout = squeeze(AN{istim}.an_sout);
	an_sout_lo = squeeze(AN{istim}.an_sout_lo);
	an_sout_hi = squeeze(AN{istim}.an_sout_hi);
	model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
end

% Analyze model output
for istim = 1:nstim
	param = stim_params{istim};
	model_output = model_outputs{istim};
	switch param.type
		case 'typMTFN'
			[~, model_MTF, ~, ~] = plotMTF(param, model_output.avIC, 0);
		case 'TIN'
			[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
		case 'SPEC_slide' % WB-TIN only for now
			[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
			%model_WB = model_WB(:,2); % testing fmincon

			% manuscript
			%[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
			%model_WB = reshape(model_WB, [1, 40]);
	end
end

%% Testing fmincon

% Normalize to match real data 
data_MTF = data(1:26);
data_WB = data(27:end);

% Minimize -1*R
r_MTF = corrcoef(model_MTF, data_MTF);
r_WB = corrcoef(model_WB, data_WB);
r_mean = (r_MTF(1,2) + r_WB(1,2))/2;
mse = -1* r_mean;

% Display
r_MTF = corrcoef(model_MTF, data_MTF);
r_WB = corrcoef(model_WB, data_WB);
fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
	x(1), x(2), x(3), x(4), x(5)) 
fprintf('R=%0.2f, RMTF=%0.2f, RWB=%0.2f\n',...
	mse, r_MTF(1,2), r_WB(1,2))

%% 

% % Match model peak rate to data peak rate
% % model_MTF = model_MTF .* (max(data_MTF)/max(model_MTF));
% % model_TIN = (model_TIN .* (max(data_TIN)/max(model_TIN)));
% % model_WB = model_WB .* (max(data_WB)/max(model_WB));
% % model = [model_MTF; model_TIN'; model_WB];
% % mse = 1/length(data) * sum((data - model).^2);
% 
% % Subtract noise alone / unmodulated in all conditions
% % data_MTF2 = data_MTF-data_MTF(1);
% % data_TIN2 = data_TIN-data_TIN(1);
% % data_WB2 = data_WB-data_WB(1);
% % model_MTF2 = model_MTF-model_MTF(1);
% % model_TIN2 = model_TIN-model_TIN(1);
% % model_WB2 = model_WB-model_WB(1);
% % model = [model_MTF2; model_TIN2'; model_WB2];
% % data = [data_MTF2; data_TIN2; data_WB2];
% % mse = 1/length(data) * sum((data - model).^2);
% 
% % Normalize to 1 and subtract noise alone
% % data_MTF2 = data_MTF./max(data_MTF)-data_MTF(1)./max(data_MTF);
% % data_TIN2 = data_TIN./max(data_TIN)-data_TIN(1)./max(data_TIN);
% % data_WB2 = data_WB./max(data_WB)-data_WB(1)./max(data_WB);
% % model_MTF2 = model_MTF./max(model_MTF)-model_MTF(1)./max(model_MTF);
% % model_TIN2 = model_TIN./max(model_TIN)-model_TIN(1)./max(model_TIN);
% % model_WB2 = model_WB./max(model_WB)-model_WB(1)./max(model_WB);
% % model = [model_MTF2; model_TIN2'; model_WB2];
% % data = [data_MTF2; data_TIN2; data_WB2];
% % mse = 1/length(data) * sum((data - model).^2);
% 
% Minimize -1*R
%model = [model_MTF; model_TIN'; model_WB];
% r_MTF = corrcoef(model_MTF, data_MTF);
% r_TIN = corrcoef(model_TIN, data_TIN);
% r_WB = corrcoef(model_WB, data_WB);
% r_mean = (r_MTF(1,2) + r_TIN(1,2) + r_WB(1,2))/3;
% mse = -1* r_mean;
%r_all = corrcoef(data, model);
%mse = -1 * r_all(1,2); 

% % Minimizing log liklihood
% % model = [model_MTF; model_TIN'; model_WB];
% % residuals = data - model;
% % n = length(data);
% % sigma = std(residuals);
% % mse = 0.5*n*log(2*pi*sigma^2) + 0.5*sum(residuals.^2)/sigma^2;
% 
% %figure; hold on; plot(data); plot(model)
% 
% % Display
% r_MTF = corrcoef(model_MTF, data_MTF);
% r_TIN = corrcoef(model_TIN, data_TIN);
% r_WB = corrcoef(model_WB, data_WB);
% fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
% 	x(1), x(2), x(3), x(4), x(5)) 
% fprintf('MSE=%0.2f, RMTF=%0.2f, RTIN=%0.2f, RWB=%0.2f\n',...
% 	mse, r_MTF(1,2), r_TIN(1,2), r_WB(1,2))

%% Manuscript
% 
% % Split data up 
% data_MTF = data(1:26);
% data_TIN = data(27:35);
% data_WB1 = data(36:75);
% data_WB2 = data(76:115);
% data_WB3 = data(116:155);
% 
% % Put model data into this format 
% model_TIN = reshape(model_TIN', [1, 9]); % not smoothed 
% for iWB = 3:5
% 	rate(iWB-2,:) = reshape(model_WB{iWB},[1, 40]);
% end
% model_WB1 = rate(1, :);
% model_WB2 = rate(2, :);
% model_WB3 = rate(3, :);
% 
% % Calculate r values 
% r_MTF = corrcoef(model_MTF, data_MTF);
% r_TIN = corrcoef(model_TIN, data_TIN);
% r_WB1 = corrcoef(model_WB1, data_WB1);
% r_WB2 = corrcoef(model_WB2, data_WB2);
% r_WB3 = corrcoef(model_WB3, data_WB3);
% r_mean = (r_MTF(1,2) + r_TIN(1,2) + r_WB1(1,2) +r_WB2(1,2)+r_WB3(1,2))/5;
% mse = -1* r_mean;
% 
% % Display
% fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
% 	x(1), x(2), x(3), x(4), x(5)) 
% fprintf('MSE=%0.2f, RMTF=%0.2f, RTIN=%0.2f, RWB1=%0.2f, RWB2=%0.2f, RWB3=%0.2f\n',...
% 	mse, r_MTF(1,2), r_TIN(1,2), r_WB1(1,2), r_WB2(1,2), r_WB3(1,2))

%% All levels 
% 
% % Split data up 
% data_MTF = data(1:26);
% data_TIN = data(27:29);
% data_WB1 = data(30:49);
% data_WB2 = data(50:69);
% data_WB3 = data(70:89);
% 
% % Calculate r values 
% r_MTF = corrcoef(model_MTF, data_MTF);
% r_TIN = corrcoef(model_TIN, data_TIN);
% r_WB1 = corrcoef(model_WB{3}(:,2), data_WB1);
% r_WB2 = corrcoef(model_WB{4}(:,2), data_WB2);
% r_WB3 = corrcoef(model_WB{5}(:,2), data_WB3);
% r_mean = (r_MTF(1,2) + r_TIN(1,2) + r_WB1(1,2) +r_WB2(1,2)+r_WB3(1,2))/5;
% mse = -1* r_mean;
% 
% % Display
% fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
% 	x(1), x(2), x(3), x(4), x(5)) 
% fprintf('MSE=%0.2f, RMTF=%0.2f, RTIN=%0.2f, RWB1=%0.2f, RWB2=%0.2f, RWB3=%0.2f\n',...
% 	mse, r_MTF(1,2), r_TIN(1,2), r_WB1(1,2), r_WB2(1,2), r_WB3(1,2))

%% One level, two SNR WBTIN conditions 

% % Split data up 
% data_MTF = data(1:26);
% data_TIN = data(27:29);
% data_WB1 = data(30:49);
% data_WB2 = data(50:69);
% 
% % Calculate r values 
% r_MTF = corrcoef(model_MTF, data_MTF);
% r_TIN = corrcoef(model_TIN, data_TIN);
% r_WB1 = corrcoef(model_WB(1:20), data_WB1);
% r_WB2 = corrcoef(model_WB(21:end), data_WB2);
% 
% r_mean = (r_MTF(1,2) + r_TIN(1,2) + r_WB1(1,2)  +r_WB2(1, 2))/4;
% mse = -1* r_mean;
% 
% % Display
% fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
% 	x(1), x(2), x(3), x(4), x(5)) 
% fprintf('MSE=%0.2f, RMTF=%0.2f, RTIN=%0.2f, RWB=%0.2f, RWB=%0.2f\n',...
% 	mse, r_MTF(1,2), r_TIN(1,2), r_WB1(1,2), r_WB2(1, 2))

end
