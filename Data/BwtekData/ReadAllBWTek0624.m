%% Read all the BWTek data files in a directory
Directory = '.\DiffuseToGlobal20160624\';
BWTekFiles = dir([Directory '*.txt']);
for iFil = 1:numel(BWTekFiles)
    if iFil == 1
      BWTekData = ReadBWTek([Directory BWTekFiles(iFil).name], ';');
    else
      BWTekData(iFil) = ReadBWTek([Directory BWTekFiles(iFil).name], ';');  
    end
end

% Sort the structure based on the datenum 
DateNums = [BWTekData.DateNum];
[DateNumsSorted, iSorted] = sort(DateNums);
BWTekDataOrdered = BWTekData(iSorted);
plot(BWTekDataOrdered(1).Wavelength,  BWTekDataOrdered(1).Reference);
save BWTekData20160624.mat BWTekDataOrdered
Wv = BWTekData(1).Wavelength;
Reference = [BWTekDataOrdered.Reference];
Dark = [BWTekDataOrdered.Dark];
plot(Wv, Reference);
plot(Wv, Dark);

%% There is a bad dark reference measurement - root it out
for iM = 1:numel(BWTekDataOrdered)
    figure;
    plot(Wv, Dark(:,iM));
    title(num2str(iM));
end
% 23, 22, 21, 20
BWTek = BWTekDataOrdered([1:19 24:end]);
Dark = [BWTek.Dark];
for iM = 1:numel(BWTek)
    figure;
    plot(Wv, Dark(:,iM));
    title(num2str(iM));
end

%% Plot some diffuse to global data
Reference = [BWTek.Reference];
Trans = [BWTek.TR1];
Light = [BWTek.Rawdata1];
Diff2Glob = (Light - Dark) ./ (Reference - Dark);

plot(Wv, Light-Dark);
plot(Wv, Reference - Dark);
plot(Wv, Diff2Glob);
for iM = 1:numel(BWTek)
    figure;
    plot(Wv, Trans(:,iM));
    title(num2str(iM));
    ylim([0 100]);
end

% Only numbers 10 and 30 seemd to have worked !!
save BWTekData20160624.mat BWTek