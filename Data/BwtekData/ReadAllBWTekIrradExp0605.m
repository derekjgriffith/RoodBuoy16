%% Read all the BWTek data files in a directory
close all
clear all
Directory = '.\IrradExp20160605\';
BWTekFiles = dir([Directory 'SP*.csv']);

iCount = 1;
for iFil = 1:numel(BWTekFiles)
    BWTekData = ReadBWTek([Directory BWTekFiles(iFil).name]);
    if isfield(BWTekData, 'IrradiancemWcm2nm1')
        if iCount == 1
          BWTekDataUnordered = BWTekData;
        else
          BWTekDataUnordered(iCount) = BWTekData;     
        end
        figure;
        plot(BWTekData.Wavelength,BWTekData.IrradiancemWcm2nm1);
        axis([350 1000 0 0.1])
        title([num2str(iCount) ' on ' BWTekData.DateStr]); 
        grid();
        
        iCount = iCount + 1;
    end
end
% Hand select measurements based on time and magnitude
iDiff = [6 8 10 12];
iGlob = [2 3 4 5 7 9 11];
iGood = [iDiff iGlob];
plot(BWTekData.Wavelength, [BWTekDataUnordered(iDiff).IrradiancemWcm2nm1]);
plot(BWTekData.Wavelength, [BWTekDataUnordered(iGlob).IrradiancemWcm2nm1]);
GlobMeanIrrad = mean([BWTekDataUnordered(iGlob).IrradiancemWcm2nm1], 2);
DiffMeanIrrad = mean([BWTekDataUnordered(iDiff).IrradiancemWcm2nm1], 2);
plot(BWTekData.Wavelength, GlobMeanIrrad, BWTekData.Wavelength, DiffMeanIrrad);
Diff2GlobRatio = DiffMeanIrrad(34:end) ./ GlobMeanIrrad(34:end); % First 33 values are nan
Wavelength = BWTekData.Wavelength(34:end);
figure;
plot(Wavelength, Diff2GlobRatio);
axis([350, 1000, 0, 0.8])
Diff2GlobRatioSmooth = fastsmooth(Diff2GlobRatio, 6, 1, 1);
figure;
plot(Wavelength, Diff2GlobRatioSmooth);
axis([350, 1000, 0, 0.8])
title('BWTek Mean Diffuse/Global Irradiance Ratio');
xlabel('Wavelength [nm]')
ylabel('Diffuse/Global Ratio')
grid;

save BWTekDataIrradExp20160605.mat BWTekDataUnordered Wavelength GlobMeanIrrad DiffMeanIrrad Diff2GlobRatio Diff2GlobRatioSmooth


