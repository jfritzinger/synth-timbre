function [base, datapath, savepath, ppi] = getPaths()

ppi = get(0, 'ScreenPixelsPerInch');
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);

% Set the 'base' filepath for creating all figures 
if ismac
	base = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre';
else
	base = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre';
end

% baseic paths for loading data and saving figures 
datapath = fullfile(base, 'data', '2025-manuscript');
savepath = fullfile(base, 'figures');

end