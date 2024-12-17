%% Figure 9: Temporal Analysis Future Directions

%% Load examples
clear

base = getPaths();
fpath = 'data/aro-2021/Figure 8';
file_names = dir(fullfile(base, fpath, 'R0*'));

for ii = 1:length(file_names)
	data(ii) = load(fullfile(base, fpath, file_names(ii).name));
end

%% Plot data

for ii = 1:length(data)

	fig = figure('Name', sprintf('R0%dS%d, TT%dN%d, DSID %d', data(ii).param.animal, ...
		data(ii).param.session, data(ii).spike_info.tetrode, data(ii).spike_info.neuron, ...
		data(ii).param.dsid));
	tiledlayout(10, 5, 'tilespacing', 'none', 'padding', 'none');

	for i_fpeak = 1:data(ii).param.stp_otc
		x = data(ii).spike_info.spike_times{i_fpeak};
		x = x(x<0.3e6 & x>0.05e6);
		y = data(ii).spike_info.spike_reps{i_fpeak};
		y = y(x<0.3e6 & x>0.05e6);

		%     subplot(3, 1, 1)
		%     plot(x/1000, y,'k.'); % spike rasters
		%     subplot(3, 1, 2)
		%     histogram(x/1000, 300);
		%     subplot(3, 1, 3)
		%     histogram(mod(x/1e6, 1/param.Delta_F), 32)

		nexttile
		histogram(mod(x/1e6*1000, 1/data(ii).param.Delta_F*1000), 32)
		%histogram(mod(x/1e6*1000, 1/data.param.fpeak_mid*1000), 32)
		ylabel(round(data(ii).param.fpeaks(i_fpeak)))
		ylim([0 100])

	end
end


