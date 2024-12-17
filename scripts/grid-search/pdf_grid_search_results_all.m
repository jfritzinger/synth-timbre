%% pdf_grid_search_results_all

% Analyze data 
params_WB = data(6:8,2); % Gets binaural WB-TIN stimuli
params_MTF = data{3,2}; % Gets binaural MTFN
data_MTF = analyzeMTF(params_MTF);
data_WB = analyzeWBTIN(params_WB, []);
data_NB = analyzeNBTIN(params_NB, CF);

% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
if ismac
	report_path = '/Volumes/WBTIN_ModelFits';
else
	report_path = 'C:\DataFiles_JBF\WBTIN_ModelFits';
end

report_name = fullfile(filepath, sprintf('%s_%d', name, CF), ...
	sprintf('%s_%d_%s.pdf', name, CF, filename));
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';
