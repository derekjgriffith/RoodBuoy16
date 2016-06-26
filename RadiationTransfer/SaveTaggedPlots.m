function Success = SaveTaggedPlots(TagDescr, ResultsFolder, BaseFilename, Revision, FileExts, TagFontProperties)
% SaveTaggedPlots : Tag a plot in the lower left corner and save
%   Detailed explanation goes here
Success = true;
ax1 = axes('Position',[0 0 1 1], 'Visible', 'off');
text(.01,0.01, TagDescr, TagFontProperties{:});
% Print using each of the requested formats
for iFileExt = 1:numel(FileExts)
  print([ResultsFolder filesep BaseFilename 'Rev' Revision '.' FileExts{iFileExt}], ['-d'  FileExts{iFileExt}]);
end

