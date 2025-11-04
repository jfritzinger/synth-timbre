%% save_model_predictions_ST.m
%
% Script that runs either SFIE single-cell, energy, or SFIE population
% models for each neuron with a response to synthetic timbre
%
% Author: J. Fritzinger
% Created: 2022-09-11; Last revision: 2024-09-16
%
% -------------------------------------------------------------------------
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);

%% Run models 

jj = 1;
for imodel = 1:2
    switch imodel
        case 1
            model_type = 'SFIE';
        case 2
            model_type = 'Lat_Inh';
    end
    %model_type = 'Energy';
    %model_type = 'SFIE_pop';


    has_data = bin200(:,1) | bin200(:,2) | bin200(:,3) | bin200(:,4);
    indices = find(has_data);
    num_index = length(indices);
    for isesh = 1:num_index
        try
            timerVal2 = tic;
        	% Error: 52

        	% Load in session
        	putative = sessions.Putative_Units{indices(isesh)};
        	CF = sessions.CF(indices(isesh));
        	load(fullfile(datapath, 'neural_data', [putative '.mat']))
        	MTF_shape = sessions.MTF{indices(isesh)};
        	params_ST = data(6:9, 2);
        	if strcmp(MTF_shape, 'BS')
        		BMF = sessions.WMF(indices(isesh));
        	elseif strcmp (MTF_shape, 'BE')
        		BMF = sessions.BMF(indices(isesh));
        	else
        		BMF = 100;
        	end

        	if strcmp(model_type, 'SFIE') || strcmp(model_type, 'SFIE_pop')
        		AN = cell(4, 1);
        		SFIE = cell(4,1);
        	elseif strcmp(model_type, 'Energy')
        		Fs = 100000;
        		gamma_param.srate = Fs;
        		energy = cell(4, 1);
        		gamma_param.fc = CF;
        	end
        	for ispl = 1:4
        		timerVal = tic;

        		% Analysis
        		if isempty(params_ST{ispl})
        			% Did not record synthetic timbre response at this level
        		else
        			data_ST = analyzeST(params_ST(ispl), CF);
        			data_ST = data_ST{1};

        			% Get spont rate
        			params_RM = data{2, 2};
        			data_RM = analyzeRM(params_RM);
        			spont = data_RM.spont;

                    % Get MTF
                    params_MTF = data{3, 2};
                    data_MTF = analyzeMTF(params_MTF);

        			% Set up stimuli
        			if strcmp(model_type, 'SFIE_pop')
        				CFs = params_ST{ispl}.fpeaks;
        				params_ST{ispl}.Fs = 100000;
        				params_ST{ispl}.physio = 1;
        				params_ST{ispl}.mnrep = 5;
        				params_ST{ispl}.dur = 0.3;
        				params_ST{ispl}.stp_otc = 1;
        				[params_ST{ispl}] = generate_ST(params_ST{ispl});
        				params_ST{ispl}.num_stim = size(params_ST{ispl}.stim, 1);
        			else
        				params_ST{ispl}.Fs = 100000;
        				params_ST{ispl}.physio = 1;
        				params_ST{ispl}.mnrep = 5;
        				params_ST{ispl}.dur = 0.3;
        				[params_ST{ispl}] = generate_ST(params_ST{ispl});
        				params_ST{ispl}.num_stim = size(params_ST{ispl}.stim, 1);
        			end

        			if strcmp(model_type, 'SFIE')
        				%%

        				% Model parameters
        				model_params.type = 'SFIE';
        				model_params.range = 2; % 1 = population model, 2 = single cell model
        				model_params.species = 1; % 1 = cat, 2 = human
        				model_params.BMF = BMF;
        				model_params.CF_range = CF;
        				model_params.num_CFs = 1;
        				model_params.CFs = CF;
        				model_params.nAN_fibers_per_CF = 10;
        				model_params.cohc = 1; % (0-1 where 1 is normal)
        				model_params.cihc = 1; % (0-1 where 1 is normal)
        				model_params.nrep = 1; % how many times to run the AN model
        				model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
        				model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
        				model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
        				model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
        				model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

        				% Run model
        				AN_temp = modelAN(params_ST{ispl}, model_params); % HSR for IC input
        				SFIE_temp = wrapperIC(AN_temp.an_sout, params_ST{ispl}, model_params); % SFIE output

        				% Run model to get spont
        				params_RM.type = 'RM';
        				params_RM.dur = 0.2;
        				params_RM.ramp_dur = 0.01;
        				params_RM.reptim = 0.6;
        				params_RM.nrep = 3;
        				params_RM.freqs = [3000,3001, 1];
        				params_RM.spls = [];
        				params_RM.binmode = 2;
        				params_RM.onsetWin = 25;
        				params_RM.mnrep = 1;
        				params_RM.Fs = 100000;
        				params_RM = generate_RM(params_RM);
        				params_RM.num_stim = size(params_RM.stim, 1);
        				AN_spont = modelAN(params_RM, model_params);
        				SFIE_spont = wrapperIC(AN_spont.an_sout, params_RM, model_params);

                        % Run MTF model
                        params_MTF.type = 'MTF';
                        params_MTF.Fs = 100000;
                        params_MTF.ramp_dur = 0.05;
                        params_MTF.noise_state = 0;
                        params_MTF.noise_band = [100, 10000];
                        params_MTF.dur = 1;
                        params_MTF.reptim = 1.5;
                        params_MTF.fms = [2, 600, 3];
                        params_MTF.mdepths = [0,0,1];
                        params_MTF.binmode = 2;
                        params_MTF.No = 30;
                        params_MTF.spl = 30;
                        params_MTF.raised_sine = 1;
                        params_MTF.onsetWin = 25;
                        params_MTF.mnrep = 3; % number of model repetitions
                        params_MTF = generate_MTF(params_MTF);
                        params_MTF.num_stim = size(params_MTF.stim, 1);
                        AN_MTF = modelAN(params_MTF, model_params);
        				SFIE_MTF = wrapperIC(AN_MTF.an_sout, params_MTF, model_params);

        				% Plot
        				if strcmp(MTF_shape, 'BS')
        					[rate, rate_std] = plotST(params_ST{ispl}, SFIE_temp.average_ic_sout_BS, 0);
        					rate = rate ./ (max(rate)/max(data_ST.rate));
        					rmse = calculateRMSE(rate, data_ST.rate);

        					SFIE{ispl}.rate = rate;
        					SFIE{ispl}.rate_std = rate_std;
        					SFIE{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        					R = corrcoef(data_ST.rate,rate); % Correlation
        					SFIE{ispl}.R = R(1, 2);
        					SFIE{ispl}.R2 = R(1, 2).^2;
        					SFIE{ispl}.PSTH = plotST_PSTH(params_ST{ispl}, SFIE_temp.ic_BS, 0);
        					SFIE{ispl}.spont = SFIE_spont.average_ic_sout_BS;
        					SFIE{ispl}.rmse = rmse;
                            SFIE{ispl}.rmse = rmse;

                            [~, avBE, stdBE, MTF_shape_calc, rate_sm] = plotMTF(params_MTF, SFIE_MTF.average_ic_sout_BS, 0);
                            SFIE{ispl}.MTF_rate = avBE;
                            SFIE{ispl}.MTF_rate_std = stdBE;
                            SFIE{ispl}.MTF_shape_calc = MTF_shape_calc;
                            SFIE{ispl}.MTF_rate_sm = rate_sm;
                            R_MTF = corrcoef(data_MTF.rate,avBE); % Correlation
        					SFIE{ispl}.R_MTF = R_MTF(1, 2);

        				elseif strcmp(MTF_shape, 'BE')
        					[rate, rate_std] = plotST(params_ST{ispl}, SFIE_temp.average_ic_sout_BE, 0);
        					rate = rate ./ (max(rate)/max(data_ST.rate));
        					rmse = calculateRMSE(rate, data_ST.rate);

        					SFIE{ispl}.rate = rate;
        					SFIE{ispl}.rate_std = rate_std;
        					SFIE{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        					R = corrcoef(data_ST.rate,rate); % Correlation
        					SFIE{ispl}.R = R(1, 2);
        					SFIE{ispl}.R2 = R(1, 2).^2;
        					SFIE{ispl}.PSTH = plotST_PSTH(params_ST{ispl}, SFIE_temp.ic_BE, 0);
        					SFIE{ispl}.spont = SFIE_spont.average_ic_sout_BE;
        					SFIE{ispl}.rmse = rmse;

                            [~, avBE, stdBE, MTF_shape_calc, rate_sm] = plotMTF(params_MTF, SFIE_MTF.average_ic_sout_BE, 0);
                            SFIE{ispl}.MTF_rate = avBE;
                            SFIE{ispl}.MTF_rate_std = stdBE;
                            SFIE{ispl}.MTF_shape_calc = MTF_shape_calc;
                            SFIE{ispl}.MTF_rate_sm = rate_sm;
                            R_MTF = corrcoef(data_MTF.rate,avBE); % Correlation
        					SFIE{ispl}.R_MTF = R_MTF(1, 2);

        				else
        					SFIE{ispl}.rate = [];
        					SFIE{ispl}.rate_std = [];
        					SFIE{ispl}.fpeaks = [];
        					SFIE{ispl}.R = [];
        					SFIE{ispl}.R2 = [];
        					SFIE{ispl}.PSTH = [];
        					SFIE{ispl}.spont = [];
        					SFIE{ispl}.rmse = [];
        				end
        				SFIE{ispl}.MTF_shape = MTF_shape;
        				SFIE{ispl}.BMF = BMF;

        				[rate, rate_std] = plotST(params_ST{ispl}, AN_temp.average_AN_sout, 0);
        				AN{ispl}.rate = rate;
        				AN{ispl}.rate_std = rate_std;
        				AN{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        				R = corrcoef(data_ST.rate,rate);
        				AN{ispl}.R = R(1, 2);
        				AN{ispl}.R2 = R(1, 2).^2;
        				AN{ispl}.PSTH = plotST_PSTH(params_ST{ispl}, AN_temp.an_sout, 0);

        			elseif strcmp(model_type, 'Energy') % Energy model
        				%%

        				stimulus = [params_ST{ispl}.stim zeros(size(params_ST{ispl}.stim,1),0.1*Fs)];
        				tvals = (1:length(stimulus))/Fs;
        				gamma_IF_reg = zeros(1,length(tvals));
        				impaired = 0; % 0 = not impaired; 1 = 'impaired'
        				pin_gamma = zeros(size(stimulus, 1), Fs*params_ST{ispl}.dur+0.1*Fs);
        				for istim = 1:size(stimulus, 1)
        					pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired, 1);
        				end
        				pin_gamma = pin_gamma(:,1:params_ST{ispl}.dur*Fs);
        				energ_out = sqrt(mean(pin_gamma.^2,2));
        				[rate, rate_std] = plotST(params_ST{ispl}, energ_out, 0);
        				R_int = corrcoef(data_ST.rate,rate);

        				% Scale 'properly'
        				max_rate = max(data_ST.rate)-spont;
        				rate = rate ./ (max(rate)/max_rate)+spont;
        				% figure
        				% hold on
        				% plot(params_ST{ispl}.fpeaks, data_ST.rate)
        				% yline(spont)
        				% plot(params_ST{ispl}.fpeaks, rate)

        				rmse = calculateRMSE(rate, data_ST.rate);

        				energy{ispl}.energ_out = energ_out;
        				energy{ispl}.rate = rate;
        				energy{ispl}.rate_std = rate_std;
        				energy{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        				energy{ispl}.R = R_int(1,2);
        				energy{ispl}.R2 =  R_int(1, 2).^2;
        				energy{ispl}.rmse = rmse;

        			elseif strcmp(model_type, 'SFIE_pop')
        				%%

        				% Model parameters
        				model_params.type = 'SFIE';
        				model_params.range = 1; % 1 = population model, 2 = single cell model
        				model_params.species = 1; % 1 = cat, 2 = human
        				model_params.BMF = BMF;
        				model_params.CFs = CFs;
        				model_params.nAN_fibers_per_CF = 10;
        				model_params.cohc = 1; % (0-1 where 1 is normal)
        				model_params.cihc = 1; % (0-1 where 1 is normal)
        				model_params.nrep = 1; % how many times to run the AN model
        				model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
        				model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
        				model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
        				model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
        				model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

        				% Run model
        				AN_temp = modelAN(params_ST{ispl}, model_params); % HSR for IC input
        				SFIE_temp = wrapperIC(AN_temp.an_sout, params_ST{ispl}, model_params); % SFIE output

        				if strcmp(MTF_shape, 'BS')
        					SFIE_pop{ispl}.rate = mean(SFIE_temp.average_ic_sout_BS, 1);
        					SFIE_pop{ispl}.rate_std = std(SFIE_temp.average_ic_sout_BS, 1);
        					SFIE_pop{ispl}.fpeaks = SFIE_temp.CFs;
        					R_int = corrcoef(data_ST.rate,SFIE_pop{ispl}.rate);
        					SFIE_pop{ispl}.R = R_int(1,2);
        					SFIE_pop{ispl}.R2 =  R_int(1, 2).^2;
        					SFIE_pop{ispl}.CFs = SFIE_temp.CFs;
        					SFIE_pop{ispl}.temporal = squeeze(mean(SFIE_temp.ic_BS, 2));

        				elseif strcmp(MTF_shape, 'BE')
        					SFIE_pop{ispl}.rate = mean(SFIE_temp.average_ic_sout_BE, 1);
        					SFIE_pop{ispl}.rate_std = std(SFIE_temp.average_ic_sout_BE, 1);
        					SFIE_pop{ispl}.fpeaks = SFIE_temp.CFs;
        					R_int = corrcoef(data_ST.rate,SFIE_pop{ispl}.rate);
        					SFIE_pop{ispl}.R = R_int(1,2);
        					SFIE_pop{ispl}.R2 =  R_int(1, 2).^2;
        					SFIE_pop{ispl}.CFs = SFIE_temp.CFs;
        					SFIE_pop{ispl}.temporal = squeeze(mean(SFIE_temp.ic_BE, 2));
        				else
        					SFIE_pop{ispl}.rate = [];
        					SFIE_pop{ispl}.rate_std = [];
        					SFIE_pop{ispl}.fpeaks = [];
        					SFIE_pop{ispl}.R = [];
        					SFIE_pop{ispl}.R2 = [];
        					SFIE_pop{ispl}.CFs =[];
        					SFIE_pop{ispl}.temporal = [];
        				end
        				SFIE_pop{ispl}.MTF_shape = MTF_shape;
        				SFIE_pop{ispl}.BMF = BMF;

        				AN_pop{ispl}.rate = mean(AN_temp.average_AN_sout, 1);
        				AN_pop{ispl}.rate_std = std(AN_temp.average_AN_sout, 1);
        				AN_pop{ispl}.fpeaks = AN_temp.CFs;
        				R_int = corrcoef(data_ST.rate,AN_pop{ispl}.rate);
        				AN_pop{ispl}.R = R_int(1,2);
        				AN_pop{ispl}.R2 =  R_int(1, 2).^2;
        				AN_pop{ispl}.CFs = AN_temp.CFs;
        				AN_pop{ispl}.temporal = squeeze(mean(AN_temp.an_sout, 1));

        			elseif strcmp(model_type, 'Lat_Inh')
        				%%

                        % 				modelpath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript\model-fits';
                        % 				foldername = putative;
                        % 				d = dir(fullfile(modelpath, foldername));
                        % 				if ~isempty(d)
                        % 					load(fullfile(modelpath, foldername, [putative '_BestModel.mat']), 'fit_params', 'AN_best')
                        % 					CS_params = [fit_params(1:2) 0.001];
                        % 					BMFs = fit_params(3:5);
                        % 					oct_range = AN_best{1}.CF_span;
                        % 				else
                        % 					continue
    					if strcmp(MTF_shape, 'BS')
    						S = 0.25; % Strength, S =
    						D = 0; % Delay, D =
    						oct_range = 0.75; % CF range =
    					else
    						S = 0.4; % Strength, S =
    						D = 0; % Delay, D =
    						oct_range = 0.5; % CF range =
    					end
    					CS_params = [S S D];
    					BMFs = [100 100 100];
                        % 				end

        				% Model parameters
        				model_params.type = 'Lateral Model';
        				model_params.range = 2; % 1 = population model, 2 = single cell model
        				model_params.species = 1; % 1 = cat, 2 = human
        				model_params.num_CFs = 1;
        				model_params.BMF = 100;
        				model_params.nAN_fibers_per_CF = 10;
        				model_params.cohc = 1; % (0-1 where 1 is normal)
        				model_params.cihc = 1; % (0-1 where 1 is normal)
        				model_params.nrep = 1; % how many times to run the AN model
        				model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
        				model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
        				model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
        				model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
        				model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
        				model_params.lateral_CF = [CF*2^(-1*oct_range), CF, CF*2^oct_range];
        				model_params.CFs = model_params.lateral_CF;
        				model_params.CF_range = model_params.CFs(2);

        				% Lateral model parameters
        				model_params.config_type = 'BS inhibited by off-CF BS';


        				% Run model
        				AN_temp = modelLateralAN(params_ST{ispl}, model_params);
                        % 				latinh_temp = modelLateralSFIE(params_ST{ispl}, model_params,...
                        % 					AN_temp.an_sout, AN_temp.an_sout_lo, AN_temp.an_sout_hi,...
                        % 					'CS_params', lm_params);
        				latinh_temp = modelLateralSFIE_BMF(params_ST{ispl}, model_params, ...
        					AN_temp.an_sout, AN_temp.an_sout_lo, AN_temp.an_sout_hi, 'CS_params', CS_params,...
        					'BMFs', BMFs);

        				% Run model to get spont
        				params_RM.type = 'RM';
        				params_RM.dur = 0.2;
        				params_RM.ramp_dur = 0.01;
        				params_RM.reptim = 0.6;
        				params_RM.nrep = 1;
        				params_RM.freqs = [3000,3001, 1];
        				params_RM.spls = [];
        				params_RM.binmode = 2;
        				params_RM.onsetWin = 25;
        				params_RM.mnrep = 1;
        				params_RM.Fs = 100000;
        				params_RM = generate_RM(params_RM);
        				params_RM.num_stim = size(params_RM.stim, 1);

        				AN_spont = modelLateralAN(params_RM, model_params);
        				SFIE_spont = modelLateralSFIE_BMF(params_RM, model_params, ...
        					AN_spont.an_sout, AN_spont.an_sout_lo, AN_spont.an_sout_hi, 'CS_params', CS_params,...
        					'BMFs', BMFs);

                        % Run MTF model
                        params_MTF.type = 'MTF';
                        params_MTF.Fs = 100000;
                        params_MTF.ramp_dur = 0.05;
                        params_MTF.noise_state = 0;
                        params_MTF.noise_band = [100, 10000];
                        params_MTF.dur = 1;
                        params_MTF.reptim = 1.5;
                        params_MTF.fms = [2, 600, 3];
                        params_MTF.mdepths = [0,0,1];
                        params_MTF.binmode = 2;
                        params_MTF.No = 30;
                        params_MTF.spl = 30;
                        params_MTF.raised_sine = 1;
                        params_MTF.onsetWin = 25;
                        params_MTF.mnrep = 3; % number of model repetitions
                        params_MTF = generate_MTF(params_MTF);
                        params_MTF.num_stim = size(params_MTF.stim, 1);
                        AN_MTF = modelLateralAN(params_MTF, model_params);
        				SFIE_MTF = modelLateralSFIE_BMF(params_MTF, model_params, ...
        					AN_MTF.an_sout, AN_MTF.an_sout_lo, AN_MTF.an_sout_hi, 'CS_params', CS_params,...
        					'BMFs', BMFs);

                        % 				figure
                        % 				tiledlayout(3, 1)
                        % 				nexttile
                        % 				plot(squeeze(AN_spont.an_sout))
                        % 				ylim([0 170])
                        % 				title('AN')
                        % 				nexttile
                        % 				plot(squeeze(SFIE_spont.ic_BE))
                        % 				title('BE')
                        % 				nexttile
                        % 				plot(squeeze(SFIE_spont.ic_BS))
                        % 				title('BS')

        				% Plot
         				if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE') %% || ~isempty(d)
        					[rate, rate_std] = plotST(params_ST{ispl}, latinh_temp.avIC, 0);
        					rate = rate ./ (max(rate)/max(data_ST.rate));
        					rmse = calculateRMSE(rate, data_ST.rate);

        					lat_inh{ispl}.rate = rate;
        					lat_inh{ispl}.rate_std = rate_std;
        					lat_inh{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        					R = corrcoef(data_ST.rate,rate); % Correlation
        					lat_inh{ispl}.R = R(1, 2);
        					lat_inh{ispl}.R2 = R(1, 2).^2;
        					lat_inh{ispl}.PSTH = plotST_PSTH(params_ST{ispl}, latinh_temp.ic, 0);
        					lat_inh{ispl}.rmse = rmse;
        					lat_inh{ispl}.spont = SFIE_spont.avIC;

                            [~, avBE, stdBE, MTF_shape_calc, rate_sm] = plotMTF(params_MTF, SFIE_MTF.avIC, 0);
                            lat_inh{ispl}.MTF_rate = avBE;
                            lat_inh{ispl}.MTF_rate_std = stdBE;
                            lat_inh{ispl}.MTF_shape_calc = MTF_shape_calc;
                            lat_inh{ispl}.MTF_rate_sm = rate_sm;
                            R_MTF = corrcoef(data_MTF.rate,avBE); % Correlation
        					lat_inh{ispl}.R_MTF = R_MTF(1, 2);

        				else
        					lat_inh{ispl}.rate = [];
        					lat_inh{ispl}.rate_std = [];
        					lat_inh{ispl}.fpeaks = [];
        					lat_inh{ispl}.R = [];
        					lat_inh{ispl}.R2 = [];
        					lat_inh{ispl}.PSTH = [];
        					lat_inh{ispl}.rmse = [];
        					lat_inh{ispl}.spont = [];
        				end
        				lat_inh{ispl}.MTF_shape = MTF_shape;
        				lat_inh{ispl}.BMF = BMF;

        				[rate, rate_std] = plotST(params_ST{ispl}, AN_temp.average_AN_sout, 0);
        				AN_lat_inh{ispl}.rate = rate;
        				AN_lat_inh{ispl}.rate_std = rate_std;
        				AN_lat_inh{ispl}.fpeaks = params_ST{ispl}.fpeaks;
        				R = corrcoef(data_ST.rate,rate);
        				AN_lat_inh{ispl}.R = R(1, 2);
        				AN_lat_inh{ispl}.R2 = R(1, 2).^2;
        				AN_lat_inh{ispl}.PSTH = plotST_PSTH(params_ST{ispl}, AN_temp.an_sout, 0);
        				AN_lat_inh{ispl}.spont = mean(AN_spont.average_AN_sout);

        			end

        			elapsedTime = toc(timerVal)/60;
        			disp([putative ' Model took ' num2str(elapsedTime) ' minutes'])
        		end
            end

            elapsedTime = toc(timerVal2)/60;
        	disp([putative ' Model (all levels) took ' num2str(elapsedTime) ' minutes'])

        	% Save model
        	filename = [putative '_' model_type '.mat'];
        	savepath = 'C:\DataFiles_JBF\Synth-Timbre\data\manuscript';
        	if strcmp(model_type, 'Energy')
        		save(fullfile(savepath, 'energy_model', filename), 'params_ST', 'energy')
        	elseif strcmp(model_type, 'SFIE')
        		save(fullfile(savepath, 'SFIE_model', filename), 'params_ST', 'AN', 'SFIE', 'model_params')
        	elseif strcmp(model_type, 'SFIE_pop')
        		save(fullfile(savepath, 'SFIE_pop_model', filename), 'params_ST', 'AN_pop', 'SFIE_pop', 'model_params')
        	elseif strcmp(model_type, 'Lat_Inh')
        		%if ~isempty(d) % TEMPORARY
        		save(fullfile(savepath, 'lat_inh_model', filename), 'params_ST', 'lat_inh', 'AN_lat_inh', 'model_params')
        		%end
        	end

        catch
            skipped_indices{jj} = sprintf('Model: %s, index = %d, Putative: %s', ...
                model_type, isesh, sessions.Putative_Units{indices(isesh)});
            jj = jj+1;
        end
    end
end

for ii = 1:jj
    disp(skipped_indices{ii})
end

%% Functions

function rmse = calculateRMSE(predicted, actual)
rmse = sqrt(mean((predicted - actual).^2));
end