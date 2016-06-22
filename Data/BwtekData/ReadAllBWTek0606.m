%% Read all the BWTek data files in a directory
Directory = '.\Irrad20160606\';
BWTekFiles = dir([Directory '*.csv']);
for iFil = 1:numel(BWTekFiles)
    if iFil == 1
      BWTekData = ReadBWTek([Directory BWTekFiles(iFil).name]);
    else
      BWTekData(iFil) = ReadBWTek([Directory BWTekFiles(iFil).name]);  
    end
end

% Sort the structure based on the datenum 
DateNums = [BWTekData.DateNum];
[DateNumsSorted, iSorted] = sort(DateNums);
BWTekDataOrdered = BWTekData(iSorted);
plot(BWTekDataOrdered(1).Wavelength,  BWTekDataOrdered(1).IrradiancemWcm2nm1);
save BWTekData20160606.mat BWTekDataOrdered


