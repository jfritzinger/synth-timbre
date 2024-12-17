clear 
%% Load in old 

filename_old = 'TimbreSessions.xlsx';
old = readtable(filename_old, 'VariableNamingRule','preserve');
old(old.Rabbit == 0,:) = [];
num_sesh = size(old,1);

%% Load in new 

filename_new = 'TimbreSessions_All2.xlsx';
new = readtable(filename_new, 'VariableNamingRule','preserve');
new(new.Rabbit == 0,:) = [];
num_new = size(new,1);

%% 

for ii = 1:num_sesh-1

	rabbit_old = old.Rabbit(ii);
	session_old = old.Session(ii);
	tt_old = old.TT(ii);
	neuron_old = old.N(ii);

	index_new = find(rabbit_old==new.Rabbit & session_old==new.Session & ...
		tt_old==new.TT & neuron_old==new.N);

	new(index_new, 5:7) = old(ii, 5:7);

end

%% 

putative = cellfun(@(n) contains(n, 'R'), new.Putative_Units);
num_not = sum(putative);
sessions_new = new(putative,:);

%% Save

writetable(sessions_new, 'Fig0_TimbreSessions2.xlsx')
