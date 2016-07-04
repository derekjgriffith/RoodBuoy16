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
OverpassDateVec = [2016 06 05 7 50 30];  % UTC
OverpassDate = datestr(OverpassDateVec, 'YYYYmmdd');
OverpassDateNum = datenum(OverpassDateVec);
ScriptName = mfilename;
Rev = ScriptName(end-2:end);  % Obtain the revision from the filename.
ResultsFolder = ['ResultsS3on' OverpassDate 'Rev' Rev];
if ~exist(ResultsFolder, 'dir')
    mkdir(ResultsFolder);
end
% # SNAP pin export table
% #
% # Product:	S3A_OL_1_EFR____20160624T074916_20160624T075216_20160625T143115_0179_005_320_3420_LN1_O_NT_001
% # Created on:	Mon Jul 04 11:36:36 GMT+02:00 2016
% 
% # Wavelength:								0.0	0.0	0.0	0.0
% Name	X	    Y	    Lon	                Lat	                Color	                        Label		 OAA	     OZA	     SAA	     SZA
% pin_1	3625.5	919.5	28.373141064206788	-25.624095398282822	java.awt.Color[r=0,g=0,b=255]	Pin 1		-82.19312	0.6830949	37.41792	59.466022

OAA = -82.19312;  % deg. Observation azimuth angle (presumably relative to north through east, satellite from dam)
OZA = 0.6830949;  % deg. Observation zenith angle (satellite zenith angle as seen from the dam)
SAA = 37.41792;  % deg. Solar azimuth angle (presumably relative to north through east)
SZA = 59.466022;  % deg. Solar zenith angle

% Note : the observation zenith angle is also the angle at which the
% observation ray strikes the water. It is therefore necessary to compute
% the sky radiance looking up from BOA at the complementary azimuth 
% which is OAA - 180 (or OAA + 180 whichever is in the range 0 to 360).
% In MODTRAN, the solar azimuth is given relative to the
% satellite azimuth, measured positive north through east.
%% Set up major adjustables parameters for MODTRAN atmosphere
% Set up user defined solar spectrum
LSUNFL = 'T';
USRSUN = 'SUNnmCEOSThuillier2005.dat';  % CEOS-endorsed solar spectrum
SolarSpectrum = 'CEOS Thuillier 20015';
% LSUNFL = '4'; % Thuiller + Kurucz 1997
GNDALT = 1.225; % km ground altitude

%% First set up the key atmospheric model parameters
% The Area-Averaged (AA) surface reflectance is computed from a set
% of S3 pixels in a radius of 1.3 km of the observation site.
% This is computed by running MODTRAN with spectrally flat surface
% reflectance at zero, 50% and 100% albedo. The oxygen absorption and
% water absorption points are not used.
% Set up basic MODTRAN satellite observation case
% MicroTOPS readings are
% SN	DATE	TIME	LATITUDE	LONGITUDE	ALTITUDE PRESSURE SZA	AM	SDCORR	TEMP	ID	SIG440	SIG500	SIG675	SIG870	SIG936	STD440	STD500	STD675	STD870	STD936	R440_500	R500_675	R675_870	R870_936
% 10572	06/24/2016	 7:49:17	-25.622	28.37	1215	896	59.56	1.968	1.033	21.5	0	392.16	455.5	734.13	534.08	384.5	0.003	0.002	0.006	0.007	0.008	0.8609	0.6205	1.3746	1.389	
% AOT440	AOT500	AOT675	AOT870	AOT936	WATER
% 0.251	0.214	0.124	0.098	0.089	0.78

AOTwv = [440 500 675 870];
AOT = [0.251	0.214	0.124	0.098];
AOT550 = interp1(AOTwv, AOT, 550, 'pchip');
% Water vapour 1.05 cm
% Surface pressure 896 mb
% Temperature 21.5 degC
% SZA 59.56
% Time : 07:49:17 UTC
CDASTM = 'b';  % Perturb boundary layer aerosol extinction
ASTMX = 0.6;
NSSALB = 4;  % Number of single scattering albedo point to read on card
AWAVLN = [0.4 0.675 0.875 1.0];
ASSALB = [0.9 0.966 0.933 0.8];  % Roughly taken from AERONET

%H2O = ????; % cm Retrieved from S3
H2O = 0.78; % From MicroTOPS
H2O = H2O * 1.0; % tweak water vapour column
H2OSTR = ['g' num2str(H2O)];
O3 = 0.2661; % atm-cm  TBR ????????????????????????????
O3 = O3 * 1.0; % Tweak O3 column
O3STR = ['a' num2str(O3)];
NSTR = 4;  % Number of streams to use for DISORT

% Set up file name of area SNAP pixels for retrieveing area-averaged 
% ground reflectance. A diameter of about 2.6 km centred on the
% irradiance station was used.
%AreaSNAPpixels = '../Data/Sentinel3/S3A_OL_1_EFR____20160605T074147_20160605T074447_20160606T174711_0180_005_049_3419_LN1_O_NT_001_geometry_MaskAARev003.txt';
%WaterSNAPpixels = '../Data/Sentinel3/S3_0605_Roodeplaat_data.txt';
RrsFile = '..\Data\Rrs\Roodeplaat_ASD_rrs.txt';

%% Rrs data unfortunately not available for this date - buoy was not operational
RrsColumns = 1 + [1 3 5 7]; % Columns measured on 2016-06-05, see ../Data/Rrs/Roodeplaat_ASD_rrs_names.txt
%% Water surface reflectance
% Set the spectral water reflectance for the specific geometry
% Water reflectance at 550 nm can be computed from the Mobley tables
% Spectral variation can be quite complex, but in the current situation
% should be quite small.
WaterReflRho = [350 0.02; 1000 0.02];

%% Now do the heavy lifting
S3toTOAusingMODTRAN;

