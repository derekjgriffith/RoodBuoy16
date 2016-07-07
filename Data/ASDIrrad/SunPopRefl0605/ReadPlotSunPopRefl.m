%% Read and plot ASD sun-lollipop experiment data in reflectance mode
% The purpose of this experiment was to measure the diffuse component of 
% horizontal spectral irradiance as a fraction of the total horizontal 
% spectral irradiance.
%
% This experiment was conducted using a pole with an obscuration on the end
% ("lollipop") to obscure the direct component of solar irradiance. The
% reference measurement (White reference) is taken with the lollipop just
% out of the way (not obscuring the sensor head) and the "reflectance"
% measurement with the lollipop shadow falling on the head.
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
print('IrradDiffuseGlobalRatio20160605.pdf', '-dpdf');
print('IrradDiffuseGlobalRatio20160605.png', '-dpng');
% Nearest MicroTOPS measurement
% AOD
% AOT440 AOT500 AOT675 AOT870 AOT936
% 0.694 0.583 0.334 0.196 0.178
% Water Vapour
% 0.96

% Save the data for later comparison to MODTRAN diffuse/global irradiance
% computation.

save IrradDiffuseGlobalRatio20160605.mat IrradRatio IrradRatioMean



