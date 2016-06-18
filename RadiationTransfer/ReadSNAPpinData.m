function DataStruct = ReadSNAPpinData(filename, group_keys, group_REs, trigger_keys)
% ReadSNAPpinData : Read a pin data file exported from ESA SNAP
% 
% Reads pixel or pin data exported from ESA SNAP.
% 
% Usage :
%  >>> DataStruct = ReadSNAPpinData(filename);
%      Or
%  >>> DataStruct = ReadSNAPpinData(filename, group_keys, group_REs);
%      Or
%  >>> DataStruct = ReadSNAPpinData(filename, group_keys, trigger_keys);
%
% If group_keys is provided, any column headers containing the specified
% keys will be grouped into a matrix of the same name within the returned
% struct. For example, if group_keys = 'radiance', all data which has a
% header including teh string 'radiance' will be collected into a
% DataStruct.radiance structure field in the output./
%
% Default trigger_keys are {'Name', 'Pixel_X'}, typically the first column
% header in a SNAP pin/pixel export. These keys trigger recognition of the
% line of headers and reading of a single block of tab-delimited data,
% which can be mixed string and numeric data. If the element at the top of
% a particular column can be converted to a numeric value, then the whole
% column is converted to numeric using str2double().
%
% Example :
% >>> a = ReadSNAPpinData('S3SolarFluxData.txt', 'solar_flux', ...
%          'solar_flux_band_([0-9]+)')
% >>> str2double([a.solar_flux_toks{:}])
%
% The above example will read pin data exported from ESA SNAP and create
% a field in the output a.solar_flux containing data from all columns
% that match the regular expression 'solar_flux_band_([0-9]+)'. A further
% field called a.solar_flux_toks will contain tokens matching the group
% ([0-9]+), which is one or more digits. These can be extracted in numeric
% form using the second statement in the example.
% 

% Default of trigger_keys
if ~exist('trigger_keys', 'var')
    trigger_keys = {'Name', 'Pixel-X'};
end
% Allow a single string instead of a cell string
if exist('group_keys', 'var')
    if ~iscellstr(group_keys)
        group_keys = cellstr(group_keys);
    end
else
    group_keys = {};
end
if ~exist('group_REs', 'var')
    group_REs = group_keys;
elseif ~iscellstr(group_REs)
    group_REs = cellstr(group_REs);
end
% Check file exists
if ~exist(filename, 'file')
    error('ReadSNAPpinData:FileNotFound', 'File not found.')
end
% Open and process the file
fid = fopen(filename, 'rt');
DataStruct = struct;
while ~feof(fid)
    line = fgetl(fid);
    [tok, rest] = strtok(line, char(9));
    if any(strcmp(tok, trigger_keys))
        % Found start of data table
        colnames = strsplit(line, char(9));
        formatSpec = [repmat('%s', 1, numel(colnames))];
        % Read tab-delimited text
        dataArray = textscan(fid, formatSpec, 'Delimiter', char(9));
        % Determine which are numeric columns by trying to convert the
        % first element
        for iCol = 1:numel(colnames)
            % Use str2num for test, which returns empty for non-numerics
            if ~isempty(str2num(dataArray{iCol}{1}))
                % Use str2double for column conversion
                dataArray{iCol} = str2double(dataArray{iCol});
            end
        end
        DataStruct.data = dataArray;
        DataStruct.colnames = colnames;
        % Scan for data matching group_keys, append if found
        for iKey = 1:numel(group_keys)
            group_data = [];
            group_data_cols = {};
            group_tokens = {};
            for iCol = 1:numel(colnames)
                REtokens = regexp(colnames(iCol), group_REs{iKey}, 'tokens');
                if ~isempty(REtokens{1})  % strfind(colnames{iCol}, group_keys{iKey})
                    group_data_cols = [group_data_cols, colnames{iCol}];
                    if isempty(group_data)
                        group_data = dataArray{iCol};
                        group_tokens = REtokens{1};
                    else
                        group_data = [group_data dataArray{iCol}];
                        group_tokens = [group_tokens REtokens{1}];
                    end
                end
            end
            DataStruct.(group_keys{iKey}) = group_data;
            DataStruct.([group_keys{iKey} '_cols']) = group_data_cols;
            DataStruct.([group_keys{iKey} '_toks']) = group_tokens;
        end
        break;
    end
end
fclose(fid);

