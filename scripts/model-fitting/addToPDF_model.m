function [img, images] = addToPDF_model(images, fig, title, size)
import mlreportgen.dom.*

% Set figure size, recommended
fig.PaperSize = size;
fig.PaperPosition = [0 0 size];
fig.Units = 'inches';
fig.Position(3:4) = size;

% Add the plot to the document
if ismac
	temppath = '/Users/jfritzinger/Documents/MATLAB/Temp_PDF';
else
	temppath = 'C:\Users\jfritzinger\Documents\MATLAB\Temp_PDF';
end
name = sprintf('%s.svg', title);
tempfilepath = fullfile(temppath, name);
print(fig, tempfilepath, '-dsvg');
img = Image(tempfilepath);
delete(fig) %delete plot figure window
images = [images {img}];

end
