function [Name,X,Y,Lon,Lat,Color1,Label,Desc,Oa01_radiance,Oa02_radiance,Oa03_radiance,Oa04_radiance,Oa05_radiance,Oa06_radiance,Oa07_radiance,Oa08_radiance,Oa09_radiance,Oa10_radiance,Oa11_radiance,Oa12_radiance,Oa13_radiance,...
    Oa14_radiance,Oa15_radiance,Oa16_radiance,Oa17_radiance,Oa18_radiance,Oa19_radiance,Oa20_radiance, all_radiance] = importSNAPpins(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [NAME,X,Y,LON,LAT,COLOR1,LABEL,DESC,OA01_RADIANCE,OA02_RADIANCE,OA03_RADIANCE,OA04_RADIANCE,OA05_RADIANCE,OA06_RADIANCE,OA07_RADIANCE,OA08_RADIANCE,OA09_RADIANCE,OA10_RADIANCE,OA11_RADIANCE,OA12_RADIANCE,OA13_RADIANCE,OA14_RADIANCE,OA15_RADIANCE,OA16_RADIANCE,OA17_RADIANCE,OA18_RADIANCE,OA19_RADIANCE,OA20_RADIANCE]
%   = importSNAPpins(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [NAME,X,Y,LON,LAT,COLOR1,LABEL,DESC,OA01_RADIANCE,OA02_RADIANCE,OA03_RADIANCE,OA04_RADIANCE,OA05_RADIANCE,OA06_RADIANCE,OA07_RADIANCE,OA08_RADIANCE,OA09_RADIANCE,OA10_RADIANCE,OA11_RADIANCE,OA12_RADIANCE,OA13_RADIANCE,OA14_RADIANCE,OA15_RADIANCE,OA16_RADIANCE,OA17_RADIANCE,OA18_RADIANCE,OA19_RADIANCE,OA20_RADIANCE]
%   = importSNAPpins(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [Name,X,Y,Lon,Lat,Color1,Label,Desc,Oa01_radiance,Oa02_radiance,Oa03_radiance,Oa04_radiance,Oa05_radiance,Oa06_radiance,Oa07_radiance,Oa08_radiance,Oa09_radiance,Oa10_radiance,Oa11_radiance,Oa12_radiance,Oa13_radiance,Oa14_radiance,Oa15_radiance,Oa16_radiance,Oa17_radiance,Oa18_radiance,Oa19_radiance,Oa20_radiance]
%   = importSNAPpins('WaterDominatedPixelsRoodeplaatS3on20160605.txt',7, 13);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/06/13 15:24:38

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 7;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [1,2,3,4,5,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]);
rawCellColumns = raw(:, [6,8]);


%% Allocate imported array to column variable names
Name = cell2mat(rawNumericColumns(:, 1));
X = cell2mat(rawNumericColumns(:, 2));
Y = cell2mat(rawNumericColumns(:, 3));
Lon = cell2mat(rawNumericColumns(:, 4));
Lat = cell2mat(rawNumericColumns(:, 5));
Color1 = rawCellColumns(:, 1);
Label = cell2mat(rawNumericColumns(:, 6));
Desc = rawCellColumns(:, 2);
Oa01_radiance = cell2mat(rawNumericColumns(:, 7));
Oa02_radiance = cell2mat(rawNumericColumns(:, 8));
Oa03_radiance = cell2mat(rawNumericColumns(:, 9));
Oa04_radiance = cell2mat(rawNumericColumns(:, 10));
Oa05_radiance = cell2mat(rawNumericColumns(:, 11));
Oa06_radiance = cell2mat(rawNumericColumns(:, 12));
Oa07_radiance = cell2mat(rawNumericColumns(:, 13));
Oa08_radiance = cell2mat(rawNumericColumns(:, 14));
Oa09_radiance = cell2mat(rawNumericColumns(:, 15));
Oa10_radiance = cell2mat(rawNumericColumns(:, 16));
Oa11_radiance = cell2mat(rawNumericColumns(:, 17));
Oa12_radiance = cell2mat(rawNumericColumns(:, 18));
Oa13_radiance = cell2mat(rawNumericColumns(:, 19));
Oa14_radiance = cell2mat(rawNumericColumns(:, 20));
Oa15_radiance = cell2mat(rawNumericColumns(:, 21));
Oa16_radiance = cell2mat(rawNumericColumns(:, 22));
Oa17_radiance = cell2mat(rawNumericColumns(:, 23));
Oa18_radiance = cell2mat(rawNumericColumns(:, 24));
Oa19_radiance = cell2mat(rawNumericColumns(:, 25));
Oa20_radiance = cell2mat(rawNumericColumns(:, 26));
all_radiance = cell2mat(rawNumericColumns(:, 7:26))';