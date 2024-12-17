%% PCA 

%% Load in data 


% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find number of neurons in each category 
% Find sessions for target MTF type
MTF_target = 'BE';
isMTF = strcmp(sessions.MTF, MTF_target);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);

%% Create matrix 

spls = {'43', '63', '73', '83'};

for ispl = 2
	has_data = bin200(:,ispl);
	indices = find(has_data);
	num_index = length(indices);

	array_z = zeros(num_index,100);
	array = zeros(num_index,100);
	matrix = zeros(num_index, 41);
	CFs = sessions.CF(indices);
	CF_names = cell(num_index, 1);

	for isesh = 1:num_index
		if bin200(indices(isesh), ispl)==1

			% Load in session
			putative = sessions.Putative_Units{indices(isesh)};
			CF = sessions.CF(indices(isesh));
			load(fullfile(datapath, [putative '.mat']))
			params_ST = data(5+ispl, 2);

			%CFs(isesh) = CF;
			CF_names{isesh} = [num2str(round(CFs(isesh))) ' Hz'];

			% Analysis
			data_ST = analyzeST(params_ST);
			data_ST = data_ST{1};
			params_RM = data{2,2};
			data_RM = analyzeRM(params_RM);			
			spont = data_RM.spont;

			% General analysis
			rate = data_ST.rate;
			rate = rate - spont;
			fpeaks = data_ST.fpeaks;
			fpeaks_re_CF = log2(fpeaks/CF);
			%fpeaks_re_CF = fpeaks/CF;

			% Align by CF (approximately)
			f = linspace(-2, 2, 100);
			[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
			[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
			f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

			% Interpolate & get z-score
			r_interp = interp1(fpeaks_re_CF, rate,f_interp, 'spline');
			z_rate = zscore(r_interp);
			array_z(isesh, f_ind(1):f_ind(2)-1) = z_rate;
			array(isesh, f_ind(1):f_ind(2)-1) = r_interp;

			matrix(isesh, 1:length(rate)) = rate;
		else
			CF_names{isesh} = [num2str(round(CFs(isesh))) ' Hz'];
		end
	end
end


%% Do PCA 

matrix = matrix(:, 1:41);
%[coeff,score,latent,tsquared,explained] = pca(matrix);
[coeff,score,latent,tsquared,explained] = pca(array, 'algorithm','als');

% Plot first 10 components 
figure('Position',[9,638,1203,234])
tiledlayout(1, 5)
for ii = 1:5
	nexttile

	%plot(f, score(ii,:))
	plot(f,coeff(:,ii),'k','Marker','o')
	title(sprintf('PCA #%d, %0.2f%', ii, explained(ii)))

end

%% Plot data 

figure
plot(f, array(1:20,:))

