function eval_models 

r_all = NaN(num_paramCF, num_paramS, 3);
MTF_shapes = cell(num_paramCF, num_paramS);
for iparamCF = 1:num_paramCF

	% Load in IC model
	filename = sprintf('%s_IC_%d.mat', putative_timbre, iparamCF);
	load(fullfile(savepath, putative_timbre, filename), ...
		'params', 'model_params', 'lateral_model', 'paramS_range',  ...
		'paramCF_range')
	num_paramS = length(paramS_range);

	for iparamS = 1:num_paramS

		% Evaluate MTF
		[fms,~,fmi] = unique(double([params{1}.mlist.fm]).');
		[~,~,MTF_shape, ~] = MTFclassification(...
			lateral_model{iparamS, 1}.avIC,fms, fmi);
		MTF_shapes{iparamCF, iparamS} = MTF_shape;

		% Data MTF 
		data_MTF = analyzeMTF(params{1});

		% Evaluate WBTIN
		for iWB = 1:3
			if ~isempty(params{1+iWB})
				[SNRs,~,si] = unique([params{1+iWB}.mlist.SNR].');
				num_SNRs = length(SNRs);
				[fpeaks,~,fi] = unique([params{1+iWB}.mlist.fpeak].');
				num_fpeaks = length(fpeaks);
				rate_size = [num_fpeaks,num_SNRs];
				[avIC,~,~,~] = accumstats({fi,si},lateral_model{iparamS, iWB+1}.avIC, rate_size);
				data_WB = analyzeWBTIN(params(1+iWB), CF);
				r_WB = corrcoef(avIC(:,2),data_WB{1}.rate(:,2));
				r_all(iparamCF, iparamS, iWB) = r_WB(1,2)^2;
			end
		end
	end
end
r2 = mean(r_all, 3, 'omitnan');
[M, idx] = max(r2, [], 'all');
[best_CF, best_S] = ind2sub([num_paramCF, num_paramS], idx);

% Save results 
results.target_MTF = data_MTF.MTF_shape;
results.r2 = r2;
results.MTF_shapes = MTF_shapes;
results.best_r2 = M;
results.best_S = paramS_range(best_S);
results.best_CF = paramCF_range(best_CF);
results.best_D = 0;
results.best_MTF = MTF_shapes{best_CF, best_S};
filename = sprintf('%s_results.mat', putative_timbre);
save(fullfile(savepath, putative_timbre, filename), 'results')

end