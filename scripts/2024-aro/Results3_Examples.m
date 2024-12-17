%% Results3_Examples
% J. Fritzinger, updated 1/16/24
clear

%% Colors

RM_colors = {'#20116B', '#5E50A9', '#A49BD0'};
bin_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
contra_colors = {'#5BAEE8', '#1267B9', '#13357D', '#060650'};
model_colors = {'#E29A63', '#BD5D16', '#984102', '#4E2201'};
font_size = 16;

%% Load in datafiles

base = getPaths();
fpath = 'data/aro-2024';

%% Load in spreadsheet and find all (BE/BS) sessions with contra & binaural




%% Plot

% Figure Properties
figure('position', [138,167,1605,788])
tiledlayout(3,4,'TileSpacing','Compact', 'TileIndexing','columnmajor', ...
	'Padding','compact');

for ii = 1:4

	% figure('position', [138,167,450,788])
	% tiledlayout(3,1,'TileSpacing','Compact', 'TileIndexing','columnmajor', ...
	% 	'Padding','compact');

	% Examples for poster
	switch ii
		case 1 % All levels, clear bin peakier than contra
			session = 'R024S478_TT2_N1';
			CF = 1300;
			ylimits = [0 90];
		case 2 % Yes, All levels, complicated but cool - mostly sharper except at 40dB
			session = 'R025S552_TT1_N1';
			CF = 2000;
			ylimits = [0 100];
		case 3 % Yes, Could be a good example after plotting other session too
			session = 'R024S484_TT2_N1';
			session2 = 'R024S485_TT2_N1';
			CF = 1500;
			ylimits = [0 90];


		case 4 % Yes, All levels, contra does not show peaks until 80dB, bin more sensitive
			% session = 'R025S562_TT3_N1';
			% CF = 2300;
			% ylimits = [0 80];
			session = 'R027S010_TT3_N1';
			CF = 770;
			ylimits = [0 50];


	end


	%% Plot binaural
	data_colors = {'#82BB95', '#03882F', '#034E1C'};

	filename = fullfile(base, fpath, [session '.mat']);
	load(filename, 'params', 'data', 'cluster', 'stims')

	has_sc = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
		strcmp(p.SPEC_slide_type,'Spectral_Centroid')&&...
		p.binmode == 2&&mod(p.fpeak_mid, 200)==0,params);

	% Plot Spectral Centroid
	clear label
	nexttile
	data_ST = data(has_sc);
	params_ST = params(has_sc);
	num_DSIDs = length(data_ST);

	label_ind = 1;
	for ind = 1:num_DSIDs
		DSID = params_ST{ind}.dsid;
		if params_ST{ind}.spl == 40
			color = bin_colors{1};
		elseif params_ST{ind}.spl == 60
			color = bin_colors{2};
		elseif params_ST{ind}.spl == 70
			color = bin_colors{3};
		elseif params_ST{ind}.spl == 80
			color = bin_colors{4};
		end

		% Plot
		hold on
		h = errorbar(data_ST{ind}.fpeaks,data_ST{ind}.rate,data_ST{ind}.rate_std/(sqrt(params_ST{ind}.nrep)),...
			'LineWidth', 1.5, 'Color', color);
		label(label_ind) = {[num2str(params_ST{ind}.spl+3) ' dB SPL']};
		label_ind = label_ind+1;
	end

	% Labels
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
	label(label_ind) = {'CF'};
	%legend(label, 'location', 'northeast', 'fontsize', 14, 'edgecolor', 'none')
	xlabel('Spectral Peak Frequency (Hz)')
	%xticklabels([])
	xlim([data_ST{ind}.fpeaks(1) data_ST{ind}.fpeaks(end)]);
	if ii == 1
		ylabel('Avg. Rate (sp/s)')
	end
	ylim(ylimits)
	set(gca,'FontSize',font_size)
	box on
	grid on
	%title([name ' Binaural']);

	bin_widths = cellfun(@(d) d.width, data_ST);

	start_bin_inds = [0.045, 0.43, 0.52, 0.9];
	for ind = 1:3
		annotation('textbox',[start_bin_inds(ii) 0.87-0.028*(ind-1) 0.0829 0.0869],...
			'Color',data_colors{ind},...
			'String',[num2str(params_ST{ind}.spl+3) ' dB SPL'],...
			'FontSize',17,...
			'EdgeColor','none');
	end

	%% Plot contra
	clear label

	nexttile
	label_ind = 1;
	if ii == 3 %% Need to keep the contra from last session though for 43

		has_sc = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
		strcmp(p.SPEC_slide_type,'Spectral_Centroid')&&...
		p.binmode == 1,params);

		data_ST = data(has_sc);
		params_ST = params(has_sc);
		num_DSIDs = length(data_ST);
		for ind = 1:num_DSIDs
			DSID = params_ST{ind}.dsid;
			if params_ST{ind}.spl == 40
				color = contra_colors{1};
			elseif params_ST{ind}.spl == 60
				color = contra_colors{2};
			elseif params_ST{ind}.spl == 70
				color = contra_colors{3};
			elseif params_ST{ind}.spl == 80
				color = contra_colors{4};
			end

			% Plot
			hold on
			h = errorbar(data_ST{ind}.fpeaks,data_ST{ind}.rate,data_ST{ind}.rate_std/(sqrt(params_ST{ind}.nrep)),...
				'LineWidth', 1.5, 'Color', color);
			label(label_ind) = {[num2str(params_ST{ind}.spl+3) ' dB SPL']};
			label_ind = label_ind+1;
			width_1 = data_ST{ind}.width;
		end


		filename = fullfile(base, fpath, [session2 '.mat']);
		load(filename, 'params', 'data', 'cluster', 'stims')

	end

	has_sc = cellfun(@(p) strcmp(p.type,'SPEC_slide')&&...
		strcmp(p.SPEC_slide_type,'Spectral_Centroid')&&...
		p.binmode == 1,params);

	data_ST = data(has_sc);
	params_ST = params(has_sc);
	num_DSIDs = length(data_ST);
	spls = cellfun(@(p) p.spl, params_ST);
	[~, order] = sort(spls);

	for ind = 1:num_DSIDs
		DSID = params_ST{order(ind)}.dsid;
		if params_ST{order(ind)}.spl == 40
			color = contra_colors{1};
		elseif params_ST{order(ind)}.spl == 60
			color = contra_colors{2};
		elseif params_ST{order(ind)}.spl == 70
			color = contra_colors{3};
		elseif params_ST{order(ind)}.spl == 80
			color = contra_colors{4};
		end

		% Plot
		hold on
		h = errorbar(data_ST{order(ind)}.fpeaks,data_ST{order(ind)}.rate,data_ST{order(ind)}.rate_std/(sqrt(params_ST{order(ind)}.nrep)),...
			'LineWidth', 1.5, 'Color', color);
		label(label_ind) = {[num2str(params_ST{order(ind)}.spl+3) ' dB SPL']};
		label_ind = label_ind+1;
	end

	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
	label(label_ind) = {'CF'};
	%legend(label, 'location', 'northeast', 'fontsize', 14, 'edgecolor', 'none')
	xlabel('Spectral Peak Frequency (Hz)')
	xlim([data_ST{ind}.fpeaks(1) data_ST{ind}.fpeaks(end)]);
	ylim(ylimits)
	if ii == 1
		ylabel('Avg. Rate (sp/s)')
	end
	
	set(gca,'FontSize',font_size)
	box on
	grid on
	%title('Contralateral');

	if ii == 3
		contra_widths = [width_1; cellfun(@(d) d.width, data_ST)];
		contra_levels = [43, 63, 83];
	else
		contra_widths = cellfun(@(d) d.width, data_ST);
		contra_levels = cellfun(@(p) p.spl+3, params_ST(order));
	end

	start_bin_inds = [0.045, 0.43, 0.52, 0.9];
	for ind = 1:3
		annotation('textbox',[start_bin_inds(ii) 0.55-0.028*(ind-1) 0.0829 0.0869],...
			'Color',contra_colors{ind},...
			'String',[num2str(contra_levels(ind)) ' dB SPL'],...
			'FontSize',17,...
			'EdgeColor','none');
	end
	

	%% Plot widths
	levels = [43, 63, 83];

	nexttile
	hold on
	plot(levels, bin_widths, 'Color', '#03882F','LineWidth', 2);
	plot(contra_levels, contra_widths, 'Color', '#1267B9','LineWidth', 2);
	set(gca,'FontSize',font_size)
	ylim([0 2500])
	xlim([40 86])
	xticks([43 63 83])
	box on
	grid on
	xlabel('Level (dB SPL)')
	if ii == 1
		ylabel('1/2 Height BW (Hz)')
	end
	legend('Diotic', 'Contra', 'Location','northwest', 'fontsize', 18, 'edgecolor', 'none')
end
