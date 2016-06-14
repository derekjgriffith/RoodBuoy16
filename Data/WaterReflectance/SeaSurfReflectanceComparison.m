%% Comparison of sea surface reflectance computations
% Compare Zemax reflectance from dispersion models with Mobley tables at
% http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
clear all

% The following are the geometry definitions for the S3 overpass on
% 2016-06-05
OAA = 104.01066;  % deg. Observation azimuth angle (presumably relative to north through east, satellite from dam)
OZA = 14.151356;  % deg. Observation zenith angle (satellite zenith angle as seen from the dam)
SAA = 38.719933;  % deg. Solar azimuth angle (presumably relative to north through east)
SZA = 59.316036;  % deg. Solar zenith angle

% Read Zemax spectral reflectance at 14 deg AOI (angle of incidence, corresponding to OZA)
[Wavelength,SReflect,PReflect,STransmit,PTransmit,SAbsorb,PAbsorb,SPhaseT,PPhaseT,SPhaseR,PPhaseR,DiattT,DiattR,RetardT,RetardR] ...
 = importZemaxReflWv('SeaWaterReflectanceZemax14.TXT');
 
SReflectSea = SReflect;
PReflectSea = PReflect;
ReflectSea = (SReflectSea + PReflectSea)/2.0;

[Wavelength,SReflect,PReflect,STransmit,PTransmit,SAbsorb,PAbsorb,SPhaseT,PPhaseT,SPhaseR,PPhaseR,DiattT,DiattR,RetardT,RetardR] ...
 = importZemaxReflWv('PureWaterReflectanceZemax14.TXT');

SReflectPure = SReflect;
PReflectPure = PReflect;
ReflectPure = (SReflectPure + PReflectPure)/2.0;

% Compute Mobley 2015 value
WindSpeed = 2;

MobleyRhoLinear = SeaSurfReflectance(OAA-SAA, OZA, SZA, WindSpeed, 'linear');
MobleyRhoCubic = SeaSurfReflectance(OAA-SAA, OZA, SZA, WindSpeed, 'cubic');
MobleyRhoSpline = SeaSurfReflectance(OAA-SAA, OZA, SZA, WindSpeed, 'spline');
MobleyRhoLinearAlt = SeaSurfReflectance(SAA-OAA+180, OZA, SZA, WindSpeed, 'linear');  % All observation azimuth angles

Wavelength = Wavelength * 1000;  % convert to nm
plot(Wavelength, [ReflectSea, SReflectSea, PReflectSea], ...
     Wavelength, [ReflectPure, SReflectPure, PReflectPure], ...
     550, MobleyRhoLinear, 'or', ... 
     550, MobleyRhoLinearAlt, 'ob');
 
 legend('Mean Sea','P Sea','S Sea','Mean Pure','P Pure','S Pure', ['Mobley RAA=' num2str(OAA-SAA, 4)], ['Mobley RAA=' num2str(SAA-OAA+180, 4)]);
 grid();
 title(['Water Surface Reflectance SZA=' num2str(SZA, 3) ...
       ', Wind=' num2str(WindSpeed, 4) ', OZA=' num2str(OZA, 4) ...
       ', RAA=' num2str(OAA-SAA, 4)]);
 xlabel('Wavelength [nm]');
 ylabel('Effective Reflectance');
 print(['S3WaterSurfReflectance20160605Wind' num2str(WindSpeed) '.pdf'], '-dpdf');




