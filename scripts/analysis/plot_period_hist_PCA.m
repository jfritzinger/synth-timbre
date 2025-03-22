%%
clear 

%% Load in spreadsheets

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', ...
	spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);


%% Create matrix of all images

MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);
%isMTF = contains(sessions.MTF, MTF_target);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.ST_43dB);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.ST_63dB);
bin200(:,3) = cellfun(@(s) contains(s, 'R'), sessions.ST_73dB);
bin200(:,4) = cellfun(@(s) contains(s, 'R'), sessions.ST_83dB);
bin200_MTF = bin200; % & isMTF;
has_data = bin200_MTF(:,1) | bin200_MTF(:,2) | bin200_MTF(:,3) | bin200_MTF(:,4);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Create matrix for each neuron
ispl = 2;
for isesh = 1:num_sessions
	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))
		allPutative{isesh} = putative;

		if bin200_MTF(ineuron, ispl)==1
			param_ST = data(5+ispl, 2);
			data_ST = analyzeST(param_ST, CF);
			data_ST = data_ST{1};
			temporal = analyzeST_Temporal(param_ST{1}, data_ST);
			if size(temporal.p_hist)<40
				continue
			else
				allTuningMats(:,:, isesh) = temporal.p_hist(1:40,:);
			end
		end
	end
end

% Get rid of 0s?
row_sum = sum(sum(allTuningMats, 1), 2);
ind_empty = find(row_sum == 0);
allTuningMats(:, :, ind_empty) = [];
allPutative(ind_empty) = [];


% Save data matrix for easy plotting later 

%% PCA VARIANCE OF TUNING ANALYSIS

[height, width, num_mats] = size(allTuningMats); %measure sizes of image stack
vectMats = reshape(allTuningMats, height*width, num_mats)';  % vectorize each tuning matrix
X_centered = vectMats - mean(vectMats, 1); %subtract off means
[coeff, score, latent] = pca(X_centered); %run PCA

% Visualize top 8 components
num_components_to_show = 16;
figure; tiledlayout(4,4);
for n = 1:num_components_to_show
	nexttile(n);
	pcMat = reshape(coeff(:,n), height, width);
	imagesc(pcMat); %clim([-0.5 0.8]);
	title(['PC ' num2str(n)]);
	colorbar;
	axis square;
end

% Plot variance explained
figure;
explained = 100*latent/sum(latent);
plot(cumsum(explained),'k-'); hold on;
plot(cumsum(explained),'k.'); hold off;
xlabel('Number of Components');
ylabel('Variance Explained (%)');
title('Cumulative Variance Explained');

%% Find images that are most aligned with each principal component
figure; tiledlayout(4,4);
for n = 1:num_components_to_show
	[~, max_idx] = max(abs(score(:,n)));
	nexttile(n)
	orig_image = reshape(allTuningMats(:,:,max_idx), height, width);
	% imagesc(orig_image); %clim([-4 8]);
	pcolor(orig_image)
	title(['Neuron most aligned with PC ' num2str(n)],'FontSize',8);
	colorbar;
	shading interp
	axis square;
	
	grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
	grayMap = flipud(grayMap);
	colormap(grayMap)
end