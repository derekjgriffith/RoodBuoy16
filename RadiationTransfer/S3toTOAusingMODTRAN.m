% Run the script for the particular overpass and conditions before running
% this script.
%% Set up initial run of MODTRAN
% Set parallel friendly mode (allows multiple MODTRAN cases to run in parallel
% on the same computer without a conflict.
% Set up wavelength intervals and resolutions 
MODStartWv = 385; % Start wavelength in MODTRAN
MODStopWv = 955; % Stop wavelength in MODTRAN, includes up to S3 b20
MODDV = 0.05; % Spectral interval in MODTRAN
MODFWHM = 0.1; % Run at 0.1 nm smoothing (convolution) in MODTRAN
MODWv = MODStartWv:MODDV:MODStopWv; % Anticipated MODTRAN output wavelengths
Mod5.ParallelFriendly(true);
% This case is visible/near-infared (VIS/NIR) wavelengths
S3 = Mod5;    % Get a completely empty case instance
% Set up name and short description
S3 = S3.SetCaseName(['S3Rood2' OverpassDate]); % The SetCaseName method is the only way to set the CaseName property
S3.CaseDescr = ['S3 Overpass at Roodeplaat Dam on ' OverpassDate];

% Note that if a card is required, ALL parameters on that card must be set,
% even if the parameters are not used.

% Set up Card 1 (mandatory - main radiative transport)
S3.MODTRN = 'M';     % MODTRAN band model
S3.SPEED = 'S';      % Slow algorithm
S3.BINARY = 'f';     % Output will be ASCII
S3.LYMOLC = ' ';     % Exclude 16 auxiliary trace gases
S3.MODEL = 3;        % Mid-latitude winter canned atmosphere 
S3.ITYPE = 3;        % Slant path to ground
S3.IEMSCT = 2;       % Compute path radiance, including solar scatter 
S3.IMULT = -1;       % Include multiple scatter, computed at H2 (target/pixel)
S3.M1 = 0;           % Temperature/pressure default to MODEL (Mid-latitude winter profile)
S3.M2 = 0;           % Water vapor defaults to MODEL profile
S3.M3 = 0;           % Ozone defaults to MODEL profile
S3.M4 = 0;           % Methane defaults to MODEL profile
S3.M5 = 0;           % Nitrous oxide defaults to MODEL profile
S3.M6 = 0;           % Carbon monoxide defaults to MODEL profile
S3.MDEF = 0;         % Default O2, NO, SO2, NO2, NH3, and HNO3 species profiles.
S3.I_RD2C = 0;       % Normal program operation - no user input for profiles
S3.NOPRNT = 0;       % Minimize printing to Tape6 output file
S3.TPTEMP = 0;       % Temperature at H2 - not important, only VIS/NIR
S3.SURREF = '0.0';   % Set surface reflectance across whole spectrum
% Note that the setting of Card 1 parameters can be accomplished in a single
% call to the Set method as follows :
% S3 = S3.Set('MODTRN', 'M','SPEED', 'S', 'BINARY', 'f', 'LYMOLC', ' ','MODEL', 3, 'ITYPE', 3, 'IEMSCT', 2, ...
%                 'IMULT', 0, 'M1', 0, 'M2', 0, 'M3', 0, 'M4', 0, 'M5', 0 , 'M6', 0, ...
%                 'MDEF', 0, 'I_RD2C', 0, 'NOPRNT', 1, 'TPTEMP', 0, 'SURREF', '0.5');


% Set up Card 1A (mandatory - main radiative transport continued)
S3.DIS = 't';        % Using DISORT multiple scattering algorithm
S3.DISAZM = 't';     % Therefore also not using azimuth dependence in DISORT
S3.DISALB = 'f';     % Don't calculate atmospheric correction data
S3.NSTR = NSTR;        % Number of streams for DISORT
S3.SFWHM = 0;        % Default solar irradiance data
S3.CO2MX = 380;      % CO2 mixing ratio, xxx ppm by volume
S3.H2OSTR = H2OSTR;  % Scale/set of water vapor amount
S3.O3STR = O3STR;    % Scale/set of ozone profile amount
S3.C_PROF = '0';     % No scaling of default molecular species profiles
S3.LSUNFL = LSUNFL;     % Read specified solar irradiance data
S3.LBMNAM = 'f';     % Don't read alternative band model file
S3.LFLTNM = 't';     % Must read filter file specified
S3.H2OAER = 'f';     % Don't bother to modify aerosol properties on the basis of H2OSTR
S3.CDTDIR = 'f';     % Data files are in the default location
S3.SOLCON = -1;      % Unity scaling of TOA solar irradiance, but apply seasonal correction
S3.CDASTM = ' ';     % No Angstrom law manipulations
% Not really necessary to set the Angstrom law data since it will not be
% used
% S3.ASTMC
% S3.ASTMX
% S3.ASTMO
% S3.AERRH
S3.NSSALB = 0;       % Use reference aerosol single-scattering albedo

% Deal with EO camera band filters
% Read the Sentinel 3 spectral response functions
S3Flt = Mod5.ReadFlt('Sentinel3SRF2011Cam4.flt');
% Drop the last band
S3Flt.Filters = S3Flt.Filters(1:end-1);
% Plot the filters
Mod5.PlotFlt(S3Flt);
% And attach the EO camera filters to the case
S3 = S3.AttachFlt(S3Flt); % This will automatically set FILTNM (Card 1A3)

% Set up Card 2 (mandatory - main aerosol and cloud options)
S3.APLUS = '  ';     % Don't use flexible aerosol manipulations
S3.IHAZE = 1;        % Rural aerosol model, visibility = 23 km (modified below)
S3.CNOVAM = ' ';     % Don't invoke NOVAM
S3.ISEASN = 0;       % Use default seasonal aerosol tweaking
S3.ARUSS = '   ';    % Don't use extended user-defined aerosol facility
S3.IVULCN = 0;       % Background stratospheric aerosol profile
S3.ICSTL = 1;        % Continental influence of maritime aerosols - not applicable to this case
S3.ICLD = 0;         % No clouds or rain
S3.IVSA = 0;         % Don't use Army Vertical Structure Algorithm for boundary layer aerosols
S3.VIS = -AOT550;    % Negative of the AOT at 550
S3.WSS = 0;          % Use default wind speed for named MODEL
S3.WHH = 0;          % Use default 24 hr average wind speed for named MODEL
S3.RAINRT = 0;       % Rain rate is zero (mm/hour), anyway no cloud/rain (ICLD)
S3.GNDALT = GNDALT;       % Target surface (H2) is at sea level

% Set up Card 3 (mandatory - Line of sight geometry)
% To define path (LOS) geometry in this case use PHI, H1 and H2 (combination 3c in manual)
S3.H1 = 0;           % Not used in this case - we are using a slant path to ground/space
S3.H2 = 0;           % km. Target pixel is at sea level
S3.ANGLE = 0;        % Not used in this case. (Zenith angle at H1)
S3.RANGE = 0;        % Not used in this case. Path length.
S3.BETA = 0;         % Not used in this case. Earth centre angle.
S3.RO = 0;           % Not used in this case. Radius of the Earth, will default to a reasonable value.
S3.LENN = 0;         % Not used in this case. Short path/long path switch.
S3.PHI = OZA;        % degrees. Zenith angle at H2 (pixel/target) to H1 (satellite camera)

% Set up Card 3A1 (Solar scattering geometry, required for IEMSCT = 2)
S3.IPARM = 12;       % Will specify relative solar azimuth angle and solar zenith angle below (PARM1 and PARM2)
S3.IPH = 2;          % Use Mie-generated internal database for aerosol phase functions (???????????)
S3.IDAY = OverpassDateVec(1:3); % Compute day number corresponding to day of overpass to adjust TOA solar irradiance
S3.ISOURC = 0;       % The Sun is the extraterrestrial source of scattered radiation

% Set up Card 3A2 (Solar scattering geometry, also required for IEMSCT = 2)
SunRelAz = SAA - OAA;
if SunRelAz < 0
    SunRelAz = SunRelAz + 360;
end
S3.PARM1 = SunRelAz;       % deg. The Sun azimuth is x deg (north through east positive) of LOS azimuth (H2 to H1)
S3.PARM2 = SZA;      % deg. Sun zenith angle at H2 (target/pixel).
S3.PARM3 = 0;        % Not used in this case.
S3.PARM4 = 0;        % Not used in this case.
S3.TIME = 0;         % Not used in this case.
S3.PSIPO = 0;        % Not used in this case.
S3.ANGLEM = 0;       % Not used in this case.
S3.G = 0;            % Not used in this case. (Henyey-Greenstein asymmetry parameter)
% An alternative way of setting up cards 3A1 and 3A2 (shortwave source scattering geometry)
% is to use the method SetScatGeom defined as follows:
%   MC = MC.SetScatGeom(IPARM, IDAY, ISOURC, PARM, TIME, PSIPO, ANGLEM)
% In this case, the call would be :
%   S3 = S3.SetScatGeom(12, [2009 11 2], 0, [50 20]);
% (TIME, PSIPO and ANGLEM are not used and will be set 0 by SetScatGeom)
% Note that IPH must still be set explicitly, as well as G.

% Set up Card 4 (mandatory - spectral range and resolution)
S3.V1 = MODStartWv;         % Start of spectral computation range in nm (see FLAGS(1))
S3.V2 = MODStopWv;         % End of spectral computation range in nm
S3.DV = MODDV;         % Spectral increment in nm
S3.FWHM = MODFWHM;         % Convolution filter width in nm
S3.YFLAG = ' ';      % Not going to generate .plt or .psc files
S3.XFLAG = ' ';      % Not going to generate .plt or .psc files
S3.FLAGS(1) = 'N';   % Use nanometres for spectral units (FLAGS(1)).
S3.FLAGS(4) = 'A';   % Put ALL radiance components in convolved data (tp7)
% Want to output fluxes as well
S3.FLAGS(7) = 't';  % Output fluxes
S3.MLFLX = 1;  % Only output fluxes at BOA and TOA
% An alternative way of setting up (most of) Card 4 is to use the SetSpectralRange method
% The call looks as follows:
%   MC = MC.SetSpectralRange(V1, V2, DV, FWHM, Units, ConvShape, FWHMisRelative)
% and in this case would look as follows :
%   S3 = S3.SetSpectralRange(350, 650, 0.1, 2, 'N');
% (ConvShape and FWHMisRelative are not required)
% Note that it will still be necessary to set FLAGS(4) as this is not done
% by SetSpectralRange. YFLAG, XFLAG and DLIMIT can be set using the
% SetPlot method.

% Set up Card 5 (mandatory - Repeat option)
S3.IRPT = 0;         % Stop program, only one sub-case in this run

% Now run the case (execute MODTRAN on the case)
S3 = S3.Run;

% Examine the file S3.tp6 to check the integrity of the run.
% The results are in the property fields S3.tp7, S3.sc7 and S3.chn
% S3.tp7 is the raw (unconvolved) radiance and transmittance data expressed as
% a function of wavenumber at full MODTRAN spectral resolution (lots of points).
% S3.chn contains the spectral channel (band) data for the camera filters.
% The convolved data as a function of wavelength in nm is in property S3.sc7.

% Plot the S3.sc7 (convolved, wavelength in nm) data
S3.PlotSc7({'SOLSCAT','SINGSCAT', 'GRNDRFLT','DRCTRFLT', 'TOTALRAD'});

% Plot some of the raw data for interest, single scattered path radiance, direct reflected and total radiance
S3.PlotTp7({'SINGSCAT', 'DRCTRFLT', 'TOTALRAD'});

% Plot a few of the channel outputs
S3.PlotChn({'PATH_TOTAL_SCAT_SOLAR','TOTAL_TRANSM_GRND_REFLECT'})
S3.PlotChn('SPECTRAL_RADIANCE')

%% Now run again, this time with multiple surface reflectances
Albedo = 0:0.25:1.0;
SURREF = cellstr(strsplit(num2str(Albedo)));
S3SurRef = S3;
S3SurRef = S3SurRef.SetCaseName(['S3Rood' OverpassDate 'SurRef']); % The SetCaseName method is the only way to set the CaseName property
S3SurRef.CaseDescr = ['Surface Reflectance Retrieval for S3 Overpass at Roodeplaat Dam on ' OverpassDate];
S3SurRef = S3SurRef.CreateSeries('SURREF', SURREF);
S3SurRef = S3SurRef.Run();
ChnRad = zeros(numel(S3SurRef(1).chn.SPECTRAL_RADIANCE), numel(S3SurRef));
for iCase = 1:numel(S3SurRef)
    ChnRad(:, iCase) = S3SurRef(iCase).chn.SPECTRAL_RADIANCE;
end
ChnRad = ChnRad * 10000 * 1000;  % Convert from W/sr/cm^2/nm to mW/sr/m^2/nm
plot(Albedo,ChnRad'); 

%% Obtain the mean area-averaged channel radiances from the S3 overpass
% Read pixels isolated from the S3 images for retrieving area-averaged
% surface reflectance.
S3AApixData = importSNAPpixels(AreaSNAPpixels);
CentreWavelengths = mean(S3AApixData.all_lambda0);
TargetChnRad = mean(S3AApixData.all_radiance);  % mW/sr/m^2/nm
% Exclude b21
TargetChnRad = TargetChnRad(1:end-1);  % Units are mW/sr/m^2/nm
% Interpolate the surface spectral albedo
RetrievedAlbedo = [];
for iChan = 1:20
  RetrievedAlbedo(iChan) = interp1(ChnRad(iChan, :), Albedo, TargetChnRad(iChan), 'linear');
end
% Want to exclude band 16 (O2 absorption) and band 20 (water vapour absorption)
WantedChannels = [1:15 17:19];
WantedRetrievedAlbedo = RetrievedAlbedo(WantedChannels);
WantedCentreWv = CentreWavelengths(WantedChannels);
plot(WantedCentreWv, WantedRetrievedAlbedo);
%% Set up run with retrieved AA surface reflectance on black target
% Will run an NSURF = 2 case with black pixel and area-averaged albedo
% Build a surface albedo structure for the area-averaged reflectance
Roode1AA(1).Filename = 'Roode1AA.dat';
Roode1AA(1).Header = ['Area-averaged reflectance at Roodeplaat Dam on ' OverpassDate ', retrieved from Sentinel 3'];
Roode1AA(1).title = '   1   Roode1AA';
Roode1AA(1).wv = [350; WantedCentreWv'; 1000];  % Add points at start and end to fully cover run spectral range
Roode1AA(1).refl = [WantedRetrievedAlbedo(1); WantedRetrievedAlbedo'; WantedRetrievedAlbedo(end)];
Roode1AA(2).Filename = 'Roode1AA.dat';
Roode1AA(2).Header = 'Black surface';
Roode1AA(2).title = '    2   Black';
Roode1AA(2).wv = [350; 1000];
Roode1AA(2).refl = [0; 0];
% Just set up the original case instead of copying
% The following statement does everything necessary to set up NSURF = 2 case
S3 = S3.AttachAlb(Roode1AA, 2, 1);
S3 = S3.Run;

%% Calculate and plot the downwelling spectral irradiance
figure;
Wv = S3.flx.Spectral;  % nm
% Compute total downwelling
GlobalBOAirrad = S3.flx.DownDiff(:,1) + S3.flx.DirectSol(:,1);
plot(Wv, GlobalBOAirrad);
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');

% !!! Note that irradiances are in W/cm^2/nm, while
% S3.sc7 radiances are in microwatts/sr/cm^2/nm

% Also plot diffuse to global
DiffuseToGlobalBOA = S3.flx.DownDiff(:,1) ./ GlobalBOAirrad;
plot(Wv, DiffuseToGlobalBOA);
title('Diffuse/Global Ratio, Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio');
grid();

%% Read in the BWTek data for comparison 
% Uncertainty over the clock - is there a clock in the instrument ?
load(['..\Data\BWtekData\BWTekData' OverpassDate '.mat']);
plot(Wv, fastsmooth(GlobalBOAirrad, 30), BWTekDataOrdered(1).Wavelength,  BWTekDataOrdered(1).IrradiancemWcm2nm1/1000, ...
    BWTekDataOrdered(end-10).Wavelength,  BWTekDataOrdered(end-10).IrradiancemWcm2nm1/1000);
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
grid();

%% Read in and plot the ASD data - definitely UTC
load(['..\Data\ASDIrrad\ASDIrradS3on' OverpassDate '.mat'])
plot(Wv, fastsmooth(GlobalBOAirrad, 120), ASDIrradMean.Wv, ASDIrradMean.RadData/10000); % Converting to W/cm^2/nm
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
xlim([350 1000])
legend('MODTRAN Smoothed', 'ASD', 'location', 'best');
grid();
%% Compute the water-leaving radiance
% First read in the R_rs values from Mark Matthews
RoodeRrsAll = dlmread('..\Data\Rrs\Roodeplaat_ASD_rrs.txt', '\t');
RrsWv = RoodeRrsAll(:,1);
Rrs06 = RoodeRrsAll(:, 1 + [1 3 5 7]);  % Select for June 5  ???????????
plot(RrsWv, Rrs06);
title(['R_{rs} for Roodeplaat Dam on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('R_{rs} [sr^{-1}]');
legend('P1','P2','P3','P4', 'location', 'best')
grid();

% Interpolate R_rs to flux/sc7 wavelength grid
Rrs06Interp = interp1(RrsWv, Rrs06, Wv, 'linear');
plot(Wv, Rrs06Interp);
title(['R_{rs} for Roodeplaat Dam on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('R_{rs} [sr^{-1}]');
legend('P1','P2','P3','P4', 'location', 'best')
grid();

%% Compute/plot water-leaving radiance by multuplying Rrs by the GlobBOAirrad
Lw = repmat(GlobalBOAirrad, 1, size(Rrs06Interp, 2)) .* Rrs06Interp;
% plot water-leaving radiance
plot(Wv, Lw);
title(['Water-leaving Radiance at BOA, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w [W/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best')
grid();
print([ResultsFolder filesep 'LwAtBOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'LwAtBOARev' Rev '.png'], '-dpng');

%% Compute/plot water-leaving radiance at TOA by multiplying by the path
% transmittance
LwTOA = repmat(S3.sc7.TRANS, 1, size(Lw, 2)) .* Lw;
plot(Wv, LwTOA);
title(['Water-leaving Radiance at TOA, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w at TOA [W/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best')
grid();

%% Compute the sky radiance as seen by reflection
% H1 = sensor, H2 = target
S3Sky = S3.AttachAlb(Roode1AA, 1); % NSURF = 1, area-averaged albedo
S3Sky.SetCaseName(['S3Rood' OverpassDate 'Sky']);
% Modify the path geometry to looking upwards using the ANGLE parameter
S3Sky.PHI = 0;  % Submit to ANGLE input
S3Sky.ANGLE = OZA;
SunRelAzToSky = OAA + 180 - SAA;
if SunRelAzToSky < 0
    SunRelAzToSky = SunRelAzToSky + 360;
end
if SunRelAzToSky > 360
    SunRelAzToSky = SunRelAzToSky - 360;
end
S3Sky.PARM1 = SunRelAzToSky; % Set the solar-relative viewing azimuth
S3Sky.Run;

% Interpolate the water reflectance to the wavelength grid
WaterReflRhoInterp = interp1(WaterReflRho(:,1), WaterReflRho(:,2), Wv, 'pchip');
WaterReflectedSkyRadiance = S3Sky.sc7.TOTALRAD .* WaterReflRhoInterp; % microwatts/sr/cm^2/nm

%% Plot the waterleaving radiance and water-reflected sky radiance together
% Remember to convert to common units of microwatts/sr/cm^2/nm
TotalRadianceBOA = Lw * 1e6 + repmat(WaterReflectedSkyRadiance, 1, size(Lw, 2));  % microwatts/sr/cm^2/nm
plot(Wv, Lw*1e6, Wv, WaterReflectedSkyRadiance);
title(['Water-leaving and Sky-Reflected Radiance at BOA, S3 on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('L_w and Sky-Reflected L at BOA [\muW/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'Sky-Reflected', 'location', 'best');
grid();
print([ResultsFolder filesep 'LwAndWaterReflLAtBOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'LwAndWaterReflLAtBOARev' Rev '.png'], '-dpng');

%% Compute and plot total radiance at TOA
TotalLTOA = TotalRadianceBOA .* repmat(S3.sc7.TRANS, 1, size(Lw, 2)) ...
                    + repmat(S3.sc7.TOTALRAD, 1, size(Lw, 2));  % microwatts/sr/cm^2/nm
plot(Wv, TotalLTOA);
title(['Total Radiance at TOA, S3 on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Total Radiance at TOA [\muW/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best');
grid();
print([ResultsFolder filesep 'TotalLatTOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'TotalLatTOARev' Rev '.png'], '-dpng');

%% Plot LwTOA over TotalLTOA
plot(Wv, 1e6 * LwTOA ./ TotalLTOA)
title(['L_w Over Total L at TOA, S3 on ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w at TOA / Total L')
legend('P1','P2','P3','P4', 'location', 'best');
grid();
print([ResultsFolder filesep 'LwOverTotalLatTOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'LwOverTotalLatTOARev' Rev '.png'], '-dpng');

%% Extract signals for S3 bands 1 to 20
% First have to expand SRFs onto common wavelength grid
S3FltInterp = Mod5.InterpFltOnto(S3Flt, Wv);
% Compute the integrals of the SRFs
S3FltInterpIntegral = trapz(Wv, S3FltInterp);

% Just do a loop rather that acrobatics with matrix dimensions
% Final Destination si ChanLTOA - channel radiance at TOA
ChanLTOA = zeros(size(S3FltInterp, 2), size(TotalLTOA, 2));
for iChan = 1:numel(S3FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = repmat(S3FltInterp(:,iChan), 1, size(TotalLTOA, 2));
    ChanLTOAProduct = ChanFilter .* TotalLTOA;
    ChanLTOA(iChan, :) = trapz(Wv, ChanLTOAProduct) / S3FltInterpIntegral(iChan);
end
% Convert to mW/sr/m^2/nm from microwatts/sr/cm^2/nm, a factor of 10
ChanLTOAmW = ChanLTOA * 10;

% Read some water dominated pixels from the S3 image
% [Name,X,Y,Lon,Lat,Color1,Label,Desc,Oa01_radiance,Oa02_radiance,Oa03_radiance,Oa04_radiance, ...
%     Oa05_radiance,Oa06_radiance,Oa07_radiance,Oa08_radiance,Oa09_radiance,Oa10_radiance, ...
%     Oa11_radiance,Oa12_radiance,Oa13_radiance,Oa14_radiance,Oa15_radiance,Oa16_radiance, ...
%     Oa17_radiance,Oa18_radiance,Oa19_radiance,Oa20_radiance, all_radiance] ...
%  = importSNAPpins(['..\Data\Sentinel3\WaterDominatedPixelsRoodeplaatS3on' OverpassDate '.txt']);
S3SNAPpixels = ReadSNAPpinData(WaterSNAPpixels, ...
    'all_radiance', 'Oa([0-9]+)_radiance');
S3SNAPpixels.all_radiance = S3SNAPpixels.all_radiance(:,1:20);  % Take only first 20 channels
ChanWv = [400.0	412.5	442.5	490.0	510.0	560.0	620.0	665.0	...
    673.75	681.25	708.75	753.75	761.25	764.375	767.5	778.75	...
    865.0	885.0	900.0	940.0];
% Plot the S3 TOA radiances with MODTRAN TOA radiances
plot(ChanWv, ChanLTOAmW, ChanWv, S3SNAPpixels.all_radiance', 'o');
title(['TOA Channel Radiance : S3 on ' OverpassDate ' at Roodeplaat']);
xlabel('Wavelength [nm]');
ylabel('Channel Radiance [mW/sr/m^2/nm]');
legend('MOD P1', 'MOD P2', 'MOD P3', 'MOD P4', 'S3');
grid()

%% Comparison of means
if size(S3SNAPpixels.all_radiance, 1) > 1  % Take mean over all pixels
    MeanS3Radiances = mean(S3SNAPpixels.all_radiance);
else
    MeanS3Radiances = S3SNAPpixels.all_radiance;
end
MeanChanLTOAmW = mean(ChanLTOAmW, 2);
plot(ChanWv, MeanChanLTOAmW, 'o-', ChanWv, MeanS3Radiances, 'x-');

%% Plot percentage errors
plot(ChanWv, 100*2*(MeanS3Radiances-MeanChanLTOAmW')./(MeanS3Radiances+MeanChanLTOAmW'), 'o');
title(['Percentage Error : S3 vs MODTRAN at TOA on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Error (Difference over Mean) [%]');
grid;
print([ResultsFolder filesep 'RelativeErrorMODTRANvsS3Rev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'RelativeErrorMODTRANvsS3Rev' Rev '.png'], '-dpng');

%% Compare the solar flux data
GlobalTOAsolirrad = S3.flx.DirectSol(:,2);  % W/cm^2/nm
GlobalTOAsolirrad = repmat(GlobalTOAsolirrad, 1, size(S3FltInterp, 2)); % W/cm^2/nm
ChanGlobTOAsolirrad = trapz(Wv, GlobalTOAsolirrad .* S3FltInterp) ./ S3FltInterpIntegral;
ChanGlobTOAsolirrad = ChanGlobTOAsolirrad * 1000 * 10000;  % Convert to mW/m^2/nm as for S3 data
ChanGlobTOAsolDNI = ChanGlobTOAsolirrad ./ cos(deg2rad(SZA));   % Get the direct normal irradiance
% Now read in the S3 data
S3SolarIrrad = ReadSNAPpinData('..\Data\Sentinel3\S3SolarFluxDataAtRoodeplaatOn20160605.txt', 'all_solar_flux', 'solar_flux_band_([0-9]+)');
% Unfortunately not in order, so determine the order
BandOrder = str2double([S3SolarIrrad.all_solar_flux_toks{:}]);
S3SolarFlux = mean(S3SolarIrrad.all_solar_flux);
% And reorder
NewOrder = sortrows([BandOrder' [1:numel(BandOrder)]']);
S3SolarFlux = S3SolarFlux(NewOrder(:,2));
%% Plot ratio of MODTRAN solar DNI at TOA to S3 product.
plot(ChanWv, ChanGlobTOAsolDNI./S3SolarFlux, 'o');
title(['Solar DNI at TOA : MODTRAN / S3 (Spectrum ' LSUNFL ' )']);
xlabel('Wavelength [nm]')
ylabel('Solar DNI MODTRAN / S3');
grid();

%% Comparison of means after correcting to S3 solar flux
% Percentage errors
MeanChanLTOAmWCorr = MeanChanLTOAmW' .* S3SolarFlux ./ ChanGlobTOAsolDNI;
plot(ChanWv, 100*2*(MeanS3Radiances-MeanChanLTOAmWCorr)./(MeanS3Radiances+MeanChanLTOAmWCorr), 'o');
title(['Flux-Corrected Percentage Error : S3 vs MODTRAN at TOA on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Error (Difference over Mean) [%]');
grid
print([ResultsFolder filesep 'RelErrorMODTRANfluxcorrectedVsS3Rev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'RelErrorMODTRANfluxcorrectedVsS3Rev' Rev '.png'], '-dpng');

%% Now attempt a retrieval of Lw at BOA
RadianceBOAWithoutLw = WaterReflectedSkyRadiance;  % microwatts/sr/cm^2/nm
LTOAWithoutLw = RadianceBOAWithoutLw .* S3.sc7.TRANS ...
                    + S3.sc7.TOTALRAD;  % microwatts/sr/cm^2/nm
%% Compute channel radiances for LTOA Without Lw
ChanLTOAWoLw = zeros(size(S3FltInterp, 2), size(LTOAWithoutLw, 2));
for iChan = 1:numel(S3FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = repmat(S3FltInterp(:,iChan), 1, size(LTOAWithoutLw, 2));
    ChanLTOAProductWoLw = ChanFilter .* LTOAWithoutLw;
    ChanLTOAWoLw(iChan, :) = trapz(Wv, ChanLTOAProductWoLw) / S3FltInterpIntegral(iChan);
end
% Convert to mW/sr/m^2/nm from microwatts/sr/cm^2/nm, a factor of 10
ChanLTOAmWWoLw = ChanLTOAWoLw * 10;

%% Compute and plot the percentage difference between S3 at TOA
% and MODTRAN at TOA without Lw
MeanChanLTOAmWWoLw = mean(ChanLTOAmWWoLw, 2);
plot(ChanWv, MeanChanLTOAmWWoLw, 'o-', ChanWv, MeanS3Radiances, 'x-');
% Correct for solar flux differences
MeanChanLTOAmWWoLwCorr = MeanChanLTOAmWWoLw' .* S3SolarFlux ./ ChanGlobTOAsolDNI;
RetrievedLwAtTOA = MeanS3Radiances - MeanChanLTOAmWWoLwCorr;
plot(ChanWv, RetrievedLwAtTOA);

%% Now back to BOA by dividing by the path transmittance
% First need to average the path transmittance over each of the S3
% bands
S3BandPathTransmittance = zeros(size(S3FltInterp, 2), size(LTOAWithoutLw, 2));
for iChan = 1:numel(S3FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = S3FltInterp(:,iChan);
    PathTransmittanceProduct = ChanFilter .* S3.sc7.TRANS;
    % Compute Weighted mean path transmittance
    S3BandPathTransmittance(iChan) = trapz(Wv, PathTransmittanceProduct) / S3FltInterpIntegral(iChan);
end
% Retrieve Lw at BOA
RetrievedLwAtBOA = RetrievedLwAtTOA ./ S3BandPathTransmittance'; % mW/sr/m^2/nm
%% Plot retrieved L_w at BOA
% Lw is in units of W/sr/cm^2/sr, so multiply by 1000 * 100 * 100
plot(Wv, Lw*100*100*1000, ChanWv, RetrievedLwAtBOA, 'ko-');  % in mW/sr/m^2/nm
title('S3 Retrieved L_w and Calculated L_w at BOA');
xlabel('Wavelength [nm]');
ylabel('L_w at BOA [mW/sr/m^2/nm]');
grid();
print([ResultsFolder filesep 'S3RetrievedAndMeasuredLwAtBOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'S3RetrievedAndMeasuredLwAtBOARev' Rev '.png'], '-dpng');


