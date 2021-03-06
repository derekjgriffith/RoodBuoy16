function DataStruct = ReadBWTek(filename, Separator)
% ReadBWTek : Read a basic BTek spectrometer text file
if ~exist('Separator', 'var')
    Separator = ',';
end
DataStruct = struct;
fid = fopen(filename, 'r');
% Process line by line
while ~feof(fid)
    line = fgetl(fid);
    [tok, rest] = strtok(line, Separator);
    if strcmp(tok, 'Date')
        theDate = rest(2:end);
        DataStruct.DateStr = theDate;
        DataStruct.DateVec = datevec(theDate);
        DataStruct.DateNum = datenum(theDate);
        continue;
    end
    if strcmp(tok, 'Pixel')
        % Determine number of fields
        nfields = numel(strfind(line, Separator));
        theFormat = repmat('%f ', 1, nfields);
        % Read the table of data
        DataTable = textscan(fid, theFormat, 'delimiter', Separator, 'EmptyValue', NaN);
        DataStruct.table = DataTable;
        TableFieldNames = strsplit(line(1:end-1), Separator);
        for iCol = 1:numel(DataTable)
            TableFieldName = regexprep(TableFieldNames{iCol},'[^a-zA-Z0-9]','');
            DataStruct.(genvarname(TableFieldName)) = DataTable{iCol};
        end
        break;
    end
    if ~isempty(rest)
        rest = rest(2:end);
        if ~isempty(str2num(rest))
            Value = str2num(rest);
        else
            Value = rest;
        end
        FieldName = genvarname(tok);
        DataStruct.(FieldName) = Value;
        
    end
    
end
fclose(fid);




end

