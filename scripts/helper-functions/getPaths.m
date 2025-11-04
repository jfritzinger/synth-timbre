function [base, datapath, savepath, ppi] = getPaths()

ppi = get(0, 'ScreenPixelsPerInch');
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);

% Set the 'base' filepath for creating all figures 
if ismac
	%base = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02-projects/synth-timbre';
    base = ['/Users/johannafritzinger/Library/CloudStorage/Dropbox-University' ...
        'ofMichigan/Johanna Fritzinger/02-projects/synth-timbre/'];
else
	base = 'C:\Users\jfritzinger\Box\02-projects\synth-timbre';
end

% baseic paths for loading data and saving figures 
% datapath = fullfile(base, 'data', '2025-manuscript-sub2');
datapath = fullfile(base, 'data', '2025-manuscript');
savepath = fullfile(base, 'figures', '2025-manuscript-sub2');

end