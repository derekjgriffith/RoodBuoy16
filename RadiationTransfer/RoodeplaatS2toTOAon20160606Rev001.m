%% MODTRAN method of taking Sentinel 2 measurements to TOA
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
%
% The Sentinel 2 process is a little more complicated in that the product
% is Level 2 reflectance. However, the solar fluxes do not seem to be given
% in the product and have to be calculated from the MODTRAN run. This is
% not such a bad thing in that it should help to eliminate any
% inconsistencies between MODTRAN solar flux and the solar flux used in the
% product processing.

%% Set up viewing and solar azimuth geometry
% The following are exact angles taken from S2 image data
% In MODTRAN H2 is the target pixel and H1 is the satellite
% Satellite was in the east as viewed from the dam
OverpassDateVec = [2016 06 05 7 42 31];  % UTC
OverpassDate = datestr(OverpassDateVec, 'YYYYmmdd');
OverpassDateNum = datenum(OverpassDateVec);
ScriptName = mfilename;
Rev = ScriptName(end-2:end);  % Obtain the revision from the filename.
ResultsFolder = ['ResultsS2on' OverpassDate 'Rev' Rev];
if ~exist(ResultsFolder, 'dir')
    mkdir(ResultsFolder);
end
OAA = 104.01066;  % deg. Observation azimuth angle (presumably relative to north through east, satellite from dam)
OZA = 14.151356;  % deg. Observation zenith angle (satellite zenith angle as seen from the dam)
SAA = 38.719933;  % deg. Solar azimuth angle (presumably relative to north through east)
SZA = 59.316036;  % deg. Solar zenith angle

% Note : the observation zenith angle is also the angle at which the
% observation ray strikes the water. It is therefore necessary to compute
% the sky radiance looking up from BOA at the complementary azimuth 
% which is OAA - 180 (or OAA + 180 whichever is in the range 0 to 360).
% In MODTRAN, the solar azimuth is given relative to the
% satellite azimuth, measured positive north through east.
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
% The Area-Averaged (AA) surface reflectance is computed from a set
% of S2 pixels in a radius of 1.3 km of the observation site.
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
AOT = [0.703 0.615 0.362 0.206];
AOT550 = interp1(AOTwv, AOT, 550);
% Water vapour 1.05 cm
% Surface pressure 892 mb
% Temperature 22 degC
% SZA 59.31
% Time : 07:42:31 UTC
CDASTM = ' ';  % No aerosol Angstrom perturbations.
NSSALB = 0;  % No single scattering albedo manipulation

%H2O = 1.30; % cm Retrieved from S2
H2O = 1.05; % From MicroTOPS
H2O = H2O * 1.0; % tweak water vapour column
H2OSTR = ['g' num2str(H2O)];
O3 = 0.2661; % atm-cm  TBR ????????????????????????????
O3 = O3 * 1.0; % Tweak O3 column
O3STR = ['a' num2str(O3)];
NSTR = 4;  % Number of streams to use for DISORT

% Set up file name of area SNAP pixels for retrieveing area-averaged 
% ground reflectance. A diameter of about 2.6 km centred on the
% irradiance station was used.
AreaSNAPpixels = '../Data/Sentinel2/MissingFilename.txt';
WaterSNAPpixels = '../Data/Sentinel2/MissingFilename.txt';
RrsFile = '..\Data\Rrs\Roodeplaat_ASD_rrs.txt';  % To be updated ??????
RrsColumns = 1 + [1 3 5 7]; % Columns measured on 2016-06-05, see ../Data/Rrs/Roodeplaat_ASD_rrs_names.txt
%% Water surface reflectance
% Set the spectral water reflectance for the specific geometry
% Water reflectance at 550 nm can be computed from the Mobley tables
% Spectral variation can be quite complex, but in the current situation
% should be quite small.
WaterReflRho = [350 0.02; 1000 0.02];

%% Now do the heavy lifting
S2toTOAusingMODTRAN;

