%% Read and plot ASD sun-lollipop experiment data in reflectance mode
close all
clear all
IrradRatio = ASD('IrradExp*.asd.ref.txt');
IrradRatio.Plot();

for iSay = 1:numel(IrradRatio)
    IrradRatio(iSay).Plot();
end

IrradRatio = IrradRatio(2:end);
IrradRatio.Plot();

IrradRatioMean = IrradRatio.Mean();
MeanDateTime = IrradRatioMean.DateTime;
IrradRatioMean.Plot(['Mean Diffuse/Global Irradiance Ratio at ' MeanDateTime ' UTC'], [350 1000]);
ylim([0 1]);
ylabel('Diffuse/Global Horizontal Irradiance');

% Nearest MicroTOPS measurement
% AOD
% AOT440 AOT500 AOT675 AOT870 AOT936
% 0.694 0.583 0.334 0.196 0.178
% Water Vapour
% 0.96

