%% Converts Sentinel 3 spectral response functions to MODTRAN .flt
% Origin is flat text files from the libRadtran project
% www.libradtran.org
% Beyond that, the providence of this data is currently unknown
% The MATLAB Mod5 package is required to run this script.
% https://github.com/derekjgriffith/matlab-modtran-5
clear all
Directory = './';
S3filterFiles = dir([Directory filesep 'sentinel3_olci_b*']);
S3filters.UnitsHeader = 'N';
S3filters.Units = 'nm';
S3filters.FileHeader = [S3filters.UnitsHeader ' Sentinel 3 ESA Oct 2012 Cam 4'];
for iS3filter = 1:numel(S3filterFiles)
    Filename = [Directory filesep S3filterFiles(iS3filter).name];
    S3filters.FilterHeaders{iS3filter} = Filename(end-2:end); 
    % Read the data, skipping header
    filterData = importdata(Filename, ' ', 4);
    % Bracket with zeros to stop MODTRAN issuing warnings
    fData = filterData.data(:, 1:2);
    fData = [[fData(1,1) - (fData(2,1) - fData(1,1)) , 0]; ...
              fData; ...
             [fData(end, 1) + (fData(end, 1) - fData(end-1, 1)) , 0]];
    S3filters.Filters{iS3filter} = fData;
end
% Plot
Mod5.PlotFlt(S3filters);
% Save to standard MODTRAN .flt "filter" function file.
Mod5.WriteFlt(S3filters, 'Sentinel3SRF2011Cam4.flt');