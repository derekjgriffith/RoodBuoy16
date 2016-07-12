%% MODTRAN method of taking Sentinel 3 measurements to TOA
clear all
close all
% The general procedure is:
% 1) Set up atmospheric model with MicroTOPS data.
% 2) Estimate mean area-averaged surface spectral reflectance.
% 3) Compute total downwelling spectral irradiance at BOA for
%    time of satellite overpass using MODTRAN. This can be
%    checked against our shoreline measurements, but there are 2 instruments and
%    they don't agree very well. Will probably revert to using some scaling to
%    obtain a result that is scale-consistent with measurements and RTC. 
% 4) Compute water-leaving spectral radiance Lw at BOA with your measured Rrs
%    data. 
% 5) Add in water surface-reflected spectral radiance for time of
%    satellite overpass with atmospheric model (1) to obtain total upwelling
%    radiance at BOA. Surprisingly, this gave me a bit of opposition last time
%    I tried it. MODTRAN does not have a water surface reflectance model, so I
%    will use a water BRDF model with libRadtran. Sanity check can be
%    performed using computation of sky radiance and multiplying by mean water
%    reflectance calculated from refractive index of pure water. 
% 6) Propagate total upwelling radiance at BOA to TOA using atmospheric model (1)
%    spectral transmittance. 
% 7) Add in atmospheric spectral path radiance for
%    atmospheric model (1) to obtain total spectral radiance at TOA. 
% 8) Compute band spectral radiances for comparison to real S2/S3 data.
%
% Retrieval of Lw follows the same procedure except that the Lw component
% is omitted before propagation to TOA. The shortfall in band radiance at
% TOA is then taken to be surviving Lw at TOA and back-propagated to BOA
% by dividing by the atmospheric path band transmittance.

%% Set up viewing and solar azimuth geometry

% The following are exact angles taken from S3 image data
% In MODTRAN H2 is the target pixel and H1 is the satellite
% Satellite was in the east as viewed from the dam
OverpassDateVec = [2016 06 06 7 16 36];  % UTC
OverpassDate = datestr(OverpassDateVec, 'YYYYmmdd');
OverpassDateNum = datenum(OverpassDateVec);
ScriptName = mfilename;
Rev = ScriptName(end-2:end);  % Obtain the revision from the filename.
ResultsFolder = ['ResultsS3on' OverpassDate 'Rev' Rev];
if ~exist(ResultsFolder, 'dir')
    mkdir(ResultsFolder);
end
OAA = 106.80172;  % deg. Observation azimuth angle (presumably relative to north through east, satellite from dam)
OZA = 50.579746;  % deg. Observation zenith angle (satellite zenith angle as seen from the dam)
SAA = 43.870384;  % deg. Solar azimuth angle (presumably relative to north through east)
SZA = 63.2531;  % deg. Solar zenith angle

% Note : the observation zenith angle is also the angle at which the
% observation ray strikes the water. It is therefore necessary to compute
% the sky radiance looking up from BOA at the complementary azimuth 
% which is OAA - 180 (or OAA + 180 whichever is in the range 0 to 360).
% In MODTRAN, the solar azimuth is given relative to the
% satellite azimuth, measured positive north through east.
%% Set up major adjustables parameters for MODTRAN atmosphere
% MODTRAN aerosol selection with IHAZE
% = 1 RURAL extinction, default VIS = 23 km.
% = 2 RURAL extinction, default VIS = 5 km.
% = 3 NAVY MARITIME extinction. Sets VIS based on wind speed and relative humidity.
% = 4 MARITIME extinction, default VIS = 23 km (LOWTRAN model).
% = 5 URBAN extinction, default VIS = 5 km.
% = 6 TROPOSPHERIC extinction, default VIS = 50 km.
% = 7 User-defined aerosol extinction coefficients. Triggers reading CARDs 2D, 2Dl and 2D2 for up to 4 altitude regions of user-defined extinction, absorption and asymmetry parameters. (This option is kept for backward compatibility; the ARUSS = 'USS' option affords greater flexibility in specifying user-defined aerosols).
% = 8 FOG1 (Advective Fog) extinction, 0.2 km VIS.
% = 9 FOG2 (Radiative Fog) extinction, 0.5 km VIS.
% = 10 DESERT extinction, sets visibility from wind speed (WSS).
IHAZE = 1; %  Actually will tune the aerosol model
IHAZEModel = 'Tuned';
% Exoatmospheric solar irradiance, 'f' for default (Kurucz 1997), then
% 1 : Kurucz 2005
% 2 : Chance + Kurucz 1997
% 3 : Cebula + Chance + Kurucz 1997
% 4 : Thuillier + Kurucz 1997
% 5 : Fontenla
% 6 : Sub-region renormalized Kurucz 1997
% 7 : Kurucz 1995 
LSUNFL = '4';
SolarSpectrum = 'Thuillier + Kurucz 1997';
GNDALT = 1.225; % km ground altitude

%% First set up the key atmospheric model parameters
% The Area-Averaged (AA) surface reflectance is computed from a set
% of S3 pixels in a radius of 1.3 km of the observation site.
% This is computed by running MODTRAN with spectrally flat surface
% reflectance at zero, 50% and 100% albedo. The oxygen absorption and
% water absorption points are not used.
% Set up basic MODTRAN satellite observation case
% MicroTOPS readings are
% AOT440 = 0.703
% AOT500 = 0.615
% AOT675 = 0.362
% AOT870 = 0.206
% AOT936 = 0.188
AOTwv = [440 500 675 870];
AOT = [0.728 0.644 0.395 0.237];
AOT550 = interp1(AOTwv, AOT, 550);
%AOT550 = 0.45; %%%%%%%% ?????????????????????????????????????????????
CDASTM = 'b';  % Perturb boundary layer aerosol extinction
ASTMX = 0.6;
NSSALB = 4;  % Number of single scattering albedo point to read on card
AWAVLN = [0.4 0.675 0.875 1.0];
ASSALB = [0.9 0.85 0.82 0.79];  % Roughly taken from AERONET
% Water vapour 0.99 cm
% Surface pressure 894 mb
% Temperature 20.5 degC
% SZA 63.29
WindSpeed = 0.5;  % m/s, probably less than 1

%H2O = 1.1473; % cm Retrieved from S3
H2O = 0.99; % From MicroTOPS
H2O = H2O * 1.0; % tweak water vapour column
H2OSTR = ['g' num2str(H2O)];
O3 = 0.2661; % atm-cm  TBR ????????????????????????????
O3 = O3 * 1.0; % Tweak O3 column
O3STR = ['a' num2str(O3)];
NSTR = 4;  % Number of streams to use for DISORT

% Set up file name of area SNAP pixels for retrieveing area-averaged 
% ground reflectance. A diameter of about 2.6 km centred on the
% irradiance station was used.
AreaSNAPpixels = '../Data/Sentinel3/S3A_OL_1_EFR____20160606T071536_20160606T071836_20160607T180921_0179_005_063_3419_LN1_O_NT_001_geometry_Mask.txt';
WaterSNAPpixels = '../Data/Sentinel3/S3_0606_Roodeplaat_data.txt';
RrsFile = '..\Data\Rrs\Roodeplaat_ASD_rrs.txt';
RrsColumns = 1 + [2 4 6 8]; % Columns measured on 2016-06-06, see ../Data/Rrs/Roodeplaat_ASD_rrs_names.txt
%% Water surface reflectance
% Set the spectral water reflectance for the specific geometry
% Water reflectance at 550 nm can be computed from the Mobley tables
% Spectral variation can be quite complex, but in the current situation
% should be quite small.
Theta_v = OZA;
Phi_v = OAA - SAA;
Rho550 = SeaSurfReflectance(Phi_v, Theta_v, SZA, WindSpeed);
WaterReflRho = [350 Rho550; 1000 Rho550];  % Currently assume wavelength-independence

%% Now do the heavy lifting
S3toTOAusingMODTRAN;

