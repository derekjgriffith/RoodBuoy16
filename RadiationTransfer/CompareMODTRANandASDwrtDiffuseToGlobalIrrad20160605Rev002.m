%% Compare horizontal spectral irradiance (diffuse/global ratio)
% Compare MODTRAN result to result measured with the ASD
% 
%
% The nearest MicroTOPS measurement at Roodeplaat was
% SN	DATE	TIME	LATITUDE	LONGITUDE	ALTITUDE	PRESSURE	SZA	AM	SDCORR	TEMP	ID	SIG440	SIG500	SIG675	SIG870	SIG936	STD440	STD500	STD675	STD870	STD936	R440_500	R500_675	R675_870	R870_936	AOT440	AOT500	AOT675	AOT870	AOT936	WATER
% 10572	06/05/2016	 9:44:46	-25.617	28.367	1225	893	48.48	1.506	1.031	25.2	0	250.23	306.42	578.15	486.83	363.63	0.002	0.002	0.003	0	0	0.8166	0.53	1.1876	1.3388	0.694	0.583	0.334	0.196	0.178	0.96

% 
% AOD
% AOT440 AOT500 AOT675 AOT870 AOT936
% 0.694 0.583 0.334 0.196 0.178
% Water Vapour
% 0.96

clear all
close all
%% Set up solar zenith/azimuth geometry

% Here is the date on which the measurement was made using the sun lollipop
MeasureDateVec = [2016 06 05 9 47 25];  % UTC
MeasureDate = datestr(MeasureDateVec, 'YYYYmmdd');
MeasureDateNum = datenum(MeasureDateVec);
ScriptName = mfilename;
Rev = ScriptName(end-2:end);  % Obtain the revision from the filename.
ResultsFolder = ['DiffToTotalASD' MeasureDate 'Rev' Rev];
if ~exist(ResultsFolder, 'dir')
    mkdir(ResultsFolder);
end

% Calculate the solar zenith angle
% Set up the parameters of the first site location (Roodeplaat - all
% available from MicroTOPS if correctly set up)
 SiteInfo.Lat = -25.617; % degrees latitude
 SiteInfo.Long = 28.367; % degrees longitude
 SiteInfo.Alt = 1225; % metres AMSL
 SiteInfo.Press = 893; % millibars atm pressure
 SiteInfo.Temp = 25.2; % degrees C
 SiteInfo.UTTDT = 2; % input times will be UTC
% Set up Ephemeris request
 EphRequests.DateTime = MeasureDateVec;
 EphRequests.Interval = 1; % days between ephemeris points
 EphRequests.NumInter = 1; % 3 points 1 day apart
 EphRequests.Object = 0; % The sun
Ephemeris = MoshierEphem(SiteInfo, EphRequests);
% Calculate solar zenith angle from topocentric altitude
SZA = 90 - Ephemeris.TopoAlt;
SAA = Ephemeris.TopoAz;
OAA = 0;
OZA = 0;
%% Set up major adjustables parameters for MODTRAN atmosphere
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
IHAZE = 1; %
IHAZEModel = 'Tuned';

% MicroTOPS readings are

AOTwv = [440 500 675 870];
AOT = [0.694 0.583 0.334 0.196];
AOT550 = interp1(AOTwv, AOT, 550, 'spline')+0.04; % Small offset to improve fit

CDASTM = 'b';  % Perturb boundary layer aerosol extinction
ASTMX = 0.82;
NSSALB = 4;  % Number of single scattering albedo point to read on card
AWAVLN = [0.4 0.675 0.875 1.0];
ASSALB = [0.88 0.82 0.82 0.79];  % Roughly taken from AERONET
ASSALB = ASSALB * 0.95;
H2O = 0.96; % From MicroTOPS
H2O = H2O * 1.0; % tweak water vapour column
H2OSTR = ['g' num2str(H2O)];
O3 = 0.2661; % atm-cm  TBR ????????????????????????????
O3 = O3 * 1.0; % Tweak O3 column
O3STR = ['a' num2str(O3)];
NSTR = 4;  % Number of streams to use for DISORT


ASDSunPopDiff2GlobFile = '..\Data\ASDIrrad\SunPopRefl0605\IrradDiffuseGlobalRatio20160605.mat';
BWTekSunPopDiff2GlobFile = '..\Data\BWtekData\BWTekDataIrradExp20160605.mat';
%% Heavy lifting
CompareMODTRANandASDwrtDiffuseToGlobalIrrad;
