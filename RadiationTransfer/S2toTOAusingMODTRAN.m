% Run the script for the particular overpass and conditions before running
% this script.
%% Set up initial run of MODTRAN
FileExts = {'pdf', 'png'};  % Save plots in these formats
TagFontProperties = {'FontSize', 6, 'FontAngle', 'italic'}; % Tag font properties in plots
% Get a Git description of the repo to use as the plot tag in lower left
% corner
[RetCode, GitDescr] = system('git describe --dirty --always --tags');
% Set parallel friendly mode (allows multiple MODTRAN cases to run in parallel
% on the same computer without a conflict.
% Set up wavelength intervals and resolutions 
MODStartWv = 410; % Start wavelength in MODTRAN, includes S2 B1
MODStopWv = 975; % Stop wavelength in MODTRAN, includes up to S2 B9
MODDV = 0.05; % Spectral interval in MODTRAN
MODFWHM = 0.1; % Run at 0.1 nm smoothing (convolution) in MODTRAN
MODWv = MODStartWv:MODDV:MODStopWv; % Anticipated MODTRAN output wavelengths
Mod5.ParallelFriendly(true);
% This case is visible/near-infared (VIS/NIR) wavelengths
S2 = Mod5;    % Get a completely empty case instance
% Set up name and short description
S2 = S2.SetCaseName(['S2Rood' OverpassDate]); % The SetCaseName method is the only way to set the CaseName property
S2.CaseDescr = ['S2 Overpass at Roodeplaat Dam on ' OverpassDate];

% Note that if a card is required, ALL parameters on that card must be set,
% even if the parameters are not used.

% Set up Card 1 (mandatory - main radiative transport)
S2.MODTRN = 'M';     % MODTRAN band model
S2.SPEED = 'S';      % Slow algorithm
S2.BINARY = 'f';     % Output will be ASCII
S2.LYMOLC = ' ';     % Exclude 16 auxiliary trace gases
S2.MODEL = 3;        % Mid-latitude winter canned atmosphere 
S2.ITYPE = 3;        % Slant path to ground
S2.IEMSCT = 4;       % Solar radiation only, no thermal scatter
S2.IMULT = -1;       % Include multiple scatter, computed at H2 (target/pixel)
S2.M1 = 0;           % Temperature/pressure default to MODEL (Mid-latitude winter profile)
S2.M2 = 0;           % Water vapor defaults to MODEL profile
S2.M3 = 0;           % Ozone defaults to MODEL profile
S2.M4 = 0;           % Methane defaults to MODEL profile
S2.M5 = 0;           % Nitrous oxide defaults to MODEL profile
S2.M6 = 0;           % Carbon monoxide defaults to MODEL profile
S2.MDEF = 0;         % Default O2, NO, SO2, NO2, NH3, and HNO3 species profiles.
S2.I_RD2C = 0;       % Normal program operation - no user input for profiles
S2.NOPRNT = 0;       % Minimize printing to Tape6 output file
S2.TPTEMP = 0;       % Temperature at H2 - not important, only VIS/NIR
S2.SURREF = '0.0';   % Set surface reflectance across whole spectrum
% Note that the setting of Card 1 parameters can be accomplished in a single
% call to the Set method as follows :
% S2 = S2.Set('MODTRN', 'M','SPEED', 'S', 'BINARY', 'f', 'LYMOLC', ' ','MODEL', 3, 'ITYPE', 3, 'IEMSCT', 2, ...
%                 'IMULT', 0, 'M1', 0, 'M2', 0, 'M3', 0, 'M4', 0, 'M5', 0 , 'M6', 0, ...
%                 'MDEF', 0, 'I_RD2C', 0, 'NOPRNT', 1, 'TPTEMP', 0, 'SURREF', '0.5');


% Set up Card 1A (mandatory - main radiative transport continued)
S2.DIS = 't';        % Using DISORT multiple scattering algorithm
S2.DISAZM = 't';     % Therefore also not using azimuth dependence in DISORT
S2.DISALB = 'f';     % Don't calculate atmospheric correction data
S2.NSTR = NSTR;        % Number of streams for DISORT
S2.SFWHM = 0;        % Default solar irradiance data
S2.CO2MX = 380;      % CO2 mixing ratio, xxx ppm by volume
S2.H2OSTR = H2OSTR;  % Scale/set of water vapor amount
S2.O3STR = O3STR;    % Scale/set of ozone profile amount
S2.C_PROF = '0';     % No scaling of default molecular species profiles
S2.LSUNFL = LSUNFL;     % Read specified solar irradiance data
S2.LBMNAM = 'f';     % Don't read alternative band model file
S2.LFLTNM = 't';     % Must read filter file specified
S2.H2OAER = 'f';     % Don't bother to modify aerosol properties on the basis of H2OSTR
S2.CDTDIR = 'f';     % Data files are in the default location
S2.SOLCON = -1;      % Unity scaling of TOA solar irradiance, but apply seasonal correction
S2.CDASTM = CDASTM;     % Angstrom law perturbation flag.
% Set up Angstrom law perturbations if provided
if exist('ASTMC', 'var')
    S2.ASTMC = ASTMC; 
end
if exist('ASTMX', 'var')
    S2.ASTMX = ASTMX;
end
if exist('ASTMO', 'var')
    S2.ASTMO = ASTMO;
end
% S2.AERRH
S2.NSSALB = NSSALB;       % Manipulate single scattering albedo if required

% Deal with card 1A1 - user-defined solar spectrum
if strcmpi(LSUNFL, 'T')
    S2.USRSUN = ['..' filesep 'DATA' filesep USRSUN];   
end

% Card 1B for single scattering albedo
if NSSALB > 0
    S2.AWAVLN = AWAVLN;   % Wavelengths for single scattering albedo in microns
    S2.ASSALB = ASSALB;   % Single scattering albdeo
end

% Deal with EO camera band filters
% Read the Sentinel 3 spectral response functions
S2Flt = Mod5.ReadFlt('Sentinel2VNIR20110909.flt');
% Plot the filters
Mod5.PlotFlt(S2Flt);
% And attach the EO camera filters to the case
S2 = S2.AttachFlt(S2Flt); % This will automatically set FILTNM (Card 1A3)

% Set up Card 2 (mandatory - main aerosol and cloud options)
S2.APLUS = '  ';     % Don't use flexible aerosol manipulations
S2.IHAZE = 1;        % Rural aerosol model, visibility = 23 km (modified below)
S2.CNOVAM = ' ';     % Don't invoke NOVAM
S2.ISEASN = 0;       % Use default seasonal aerosol tweaking
S2.ARUSS = '   ';    % Don't use extended user-defined aerosol facility
S2.IVULCN = 0;       % Background stratospheric aerosol profile
S2.ICSTL = 1;        % Continental influence of maritime aerosols - not applicable to this case
S2.ICLD = 0;         % No clouds or rain
S2.IVSA = 0;         % Don't use Army Vertical Structure Algorithm for boundary layer aerosols
S2.VIS = -AOT550;    % Negative of the AOT at 550
S2.WSS = 0;          % Use default wind speed for named MODEL
S2.WHH = 0;          % Use default 24 hr average wind speed for named MODEL
S2.RAINRT = 0;       % Rain rate is zero (mm/hour), anyway no cloud/rain (ICLD)
S2.GNDALT = GNDALT;       % Target surface (H2) is at sea level

% Set up Card 3 (mandatory - Line of sight geometry)
% To define path (LOS) geometry in this case use PHI, H1 and H2 (combination 3c in manual)
S2.H1 = 0;           % Not used in this case - we are using a slant path to ground/space
S2.H2 = 0;           % km. Target pixel is at sea level
S2.ANGLE = 0;        % Not used in this case. (Zenith angle at H1)
S2.RANGE = 0;        % Not used in this case. Path length.
S2.BETA = 0;         % Not used in this case. Earth centre angle.
S2.RO = 0;           % Not used in this case. Radius of the Earth, will default to a reasonable value.
S2.LENN = 0;         % Not used in this case. Short path/long path switch.
S2.PHI = OZA;        % degrees. Zenith angle at H2 (pixel/target) to H1 (satellite camera)

% Set up Card 3A1 (Solar scattering geometry, required for IEMSCT = 2)
S2.IPARM = 12;       % Will specify relative solar azimuth angle and solar zenith angle below (PARM1 and PARM2)
S2.IPH = 2;          % Use Mie-generated internal database for aerosol phase functions (???????????)
S2.IDAY = OverpassDateVec(1:3); % Compute day number corresponding to day of overpass to adjust TOA solar irradiance
S2.ISOURC = 0;       % The Sun is the extraterrestrial source of scattered radiation

% Set up Card 3A2 (Solar scattering geometry, also required for IEMSCT = 2)
SunRelAz = SAA - OAA;
if SunRelAz < 0
    SunRelAz = SunRelAz + 360;
end
S2.PARM1 = SunRelAz;       % deg. The Sun azimuth is x deg (north through east positive) of LOS azimuth (H2 to H1)
S2.PARM2 = SZA;      % deg. Sun zenith angle at H2 (target/pixel).
S2.PARM3 = 0;        % Not used in this case.
S2.PARM4 = 0;        % Not used in this case.
S2.TIME = 0;         % Not used in this case.
S2.PSIPO = 0;        % Not used in this case.
S2.ANGLEM = 0;       % Not used in this case.
S2.G = 0;            % Not used in this case. (Henyey-Greenstein asymmetry parameter)
% An alternative way of setting up cards 3A1 and 3A2 (shortwave source scattering geometry)
% is to use the method SetScatGeom defined as follows:
%   MC = MC.SetScatGeom(IPARM, IDAY, ISOURC, PARM, TIME, PSIPO, ANGLEM)
% In this case, the call would be :
%   S2 = S2.SetScatGeom(12, [2009 11 2], 0, [50 20]);
% (TIME, PSIPO and ANGLEM are not used and will be set 0 by SetScatGeom)
% Note that IPH must still be set explicitly, as well as G.

% Set up Card 4 (mandatory - spectral range and resolution)
S2.V1 = MODStartWv;         % Start of spectral computation range in nm (see FLAGS(1))
S2.V2 = MODStopWv;         % End of spectral computation range in nm
S2.DV = MODDV;         % Spectral increment in nm
S2.FWHM = MODFWHM;         % Convolution filter width in nm
S2.YFLAG = ' ';      % Not going to generate .plt or .psc files
S2.XFLAG = ' ';      % Not going to generate .plt or .psc files
S2.FLAGS(1) = 'N';   % Use nanometres for spectral units (FLAGS(1)).
S2.FLAGS(4) = 'A';   % Put ALL radiance components in convolved data (tp7)
% Want to output fluxes as well
S2.FLAGS(7) = 't';  % Output fluxes
S2.MLFLX = 1;  % Only output fluxes at BOA and TOA
% An alternative way of setting up (most of) Card 4 is to use the SetSpectralRange method
% The call looks as follows:
%   MC = MC.SetSpectralRange(V1, V2, DV, FWHM, Units, ConvShape, FWHMisRelative)
% and in this case would look as follows :
%   S2 = S2.SetSpectralRange(350, 650, 0.1, 2, 'N');
% (ConvShape and FWHMisRelative are not required)
% Note that it will still be necessary to set FLAGS(4) as this is not done
% by SetSpectralRange. YFLAG, XFLAG and DLIMIT can be set using the
% SetPlot method.

% Set up Card 5 (mandatory - Repeat option)
S2.IRPT = 0;         % Stop program, only one sub-case in this run

% Now run the case (execute MODTRAN on the case)
S2 = S2.Run;

% Examine the file S2.tp6 to check the integrity of the run.
% The results are in the property fields S2.tp7, S2.sc7 and S2.chn
% S2.tp7 is the raw (unconvolved) radiance and transmittance data expressed as
% a function of wavenumber at full MODTRAN spectral resolution (lots of points).
% S2.chn contains the spectral channel (band) data for the camera filters.
% The convolved data as a function of wavelength in nm is in property S2.sc7.

% Plot the S2.sc7 (convolved, wavelength in nm) data
S2.PlotSc7({'SOLSCAT','SINGSCAT', 'GRNDRFLT','DRCTRFLT', 'TOTALRAD'});

% Plot some of the raw data for interest, single scattered path radiance, direct reflected and total radiance
S2.PlotTp7({'SINGSCAT', 'DRCTRFLT', 'TOTALRAD'});

% Plot a few of the channel outputs
S2.PlotChn({'PATH_TOTAL_SCAT_SOLAR','TOTAL_TRANSM_GRND_REFLECT'});
S2.PlotChn('SPECTRAL_RADIANCE');

%% Run with vertical path to verify the AOD/AOT
S2Trans = S2;
S2Trans.SetCaseName(['S2RoodTrans' OverpassDate]);
% Set path vertical
S2Trans.PHI = 0;
S2Trans.IEMSCT = 0;  % Switch to spectral transmittance only mode
S2Trans = S2Trans.Run;
% Record the optical depth
S2TransTotalOD = -log(S2Trans.sc7.COMBINTRANS);
% Switch off aerosols and run again
S2Trans.IHAZE = 0;
S2Trans = S2Trans.Run;
S2TransNoAerosolOD = -log(S2Trans.sc7.COMBINTRANS);
S2AOD = S2TransTotalOD - S2TransNoAerosolOD;
AODWv = S2Trans.sc7.WAVLNM;
figure;
plot(AODWv, S2AOD, AOTwv, AOT, 'o', 550, AOT550, 'o');
title(['Vertical Aerosol Optical Depth, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('AOD')
legend('MODTRAN', 'MicroTOPS', 'MicroTOPS Interpolated', 'location', 'best')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'AOD', Rev, FileExts, TagFontProperties)

%% Now run again, this time with multiple surface reflectances
Albedo = 0:0.25:1.0;
SURREF = cellstr(strsplit(num2str(Albedo)));
S2SurRef = S2;
S2SurRef = S2SurRef.SetCaseName(['S2Rood' OverpassDate 'SurRef']); % The SetCaseName method is the only way to set the CaseName property
S2SurRef.CaseDescr = ['Surface Reflectance Retrieval for S2 Overpass at Roodeplaat Dam on ' OverpassDate];
S2SurRef = S2SurRef.CreateSeries('SURREF', SURREF);
S2SurRef = S2SurRef.Run();
ChnRad = zeros(numel(S2SurRef(1).chn.SPECTRAL_RADIANCE), numel(S2SurRef));
for iCase = 1:numel(S2SurRef)
    ChnRad(:, iCase) = S2SurRef(iCase).chn.SPECTRAL_RADIANCE;
end
ChnRad = ChnRad * 10000 * 1000;  % Convert from W/sr/cm^2/nm to mW/sr/m^2/nm
%plot(Albedo, ChnRad'); 

%% Obtain the mean area-averaged channel radiances from the S2 overpass
% Read pixels isolated from the S2 images for retrieving area-averaged
% surface reflectance.
S2AApixData = ReadSNAPpinData(AreaSNAPpixels, ...
        'all_refl', 'B([0-9]+[A]?)');
CentreWavelengths = [443 490 560 665 705 740 783 842 865 945]; % B1 .. B9
%% Convert S2 reflectances to radiances
% First have to expand SRFs onto common wavelength grid
Wv = S2.flx.Spectral;  % nm
S2FltInterp = Mod5.InterpFltOnto(S2Flt, Wv);
% Compute the integrals of the SRFs
S2FltInterpIntegral = trapz(Wv, S2FltInterp);
% Determine TOA solar band irradiance
%
GlobalTOAsolirrad = S2.flx.DirectSol(:,2);  % W/cm^2/nm
GlobalTOAsolirrad = repmat(GlobalTOAsolirrad, 1, size(S2FltInterp, 2)); % W/cm^2/nm
ChanGlobTOAsolirrad = trapz(Wv, GlobalTOAsolirrad .* S2FltInterp) ./ S2FltInterpIntegral;
ChanGlobTOAsolirrad = ChanGlobTOAsolirrad * 1000 * 10000;  % Convert to mW/m^2/nm as for S2 data
ChanGlobTOAsolDNI = ChanGlobTOAsolirrad ./ cos(deg2rad(SZA));   % Get the direct normal irradiance

TargetChnRad = ChanGlobTOAsolirrad .* mean(S2AApixData.all_refl(:, 1:10))/10000/pi;  % mW/sr/m^2/nm

% Interpolate the surface spectral albedo
RetrievedAlbedo = [];
for iChan = 1:10
  RetrievedAlbedo(iChan) = interp1(ChnRad(iChan, :), Albedo, TargetChnRad(iChan), 'linear');
end
% Want to exclude band 10 (B9) (water vapour absorption)
WantedChannels = [1:9];
WantedRetrievedAlbedo = RetrievedAlbedo(WantedChannels);
WantedCentreWv = CentreWavelengths(WantedChannels);
figure;
plot(WantedCentreWv, WantedRetrievedAlbedo);
title(['Retrieved Area-Averaged Surface Reflectance, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('Area-Averaged Albedo')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'AreaAveAlbedo', Rev, FileExts, TagFontProperties);

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
S2 = S2.AttachAlb(Roode1AA, 2, 1);
S2 = S2.Run;

%% Calculate and plot the downwelling spectral irradiance
figure;

% Compute total downwelling
GlobalBOAirrad = S2.flx.DownDiff(:,1) + S2.flx.DirectSol(:,1);
figure;
plot(Wv, GlobalBOAirrad);
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
SaveTaggedPlots(GitDescr, ResultsFolder,  'TotalDownIrradBOA', Rev, FileExts, TagFontProperties);
% !!! Note that irradiances are in W/cm^2/nm, while
% S2.sc7 radiances are in microwatts/sr/cm^2/nm

%% Also plot diffuse to global
DiffuseToGlobalBOA = S2.flx.DownDiff(:,1) ./ GlobalBOAirrad;
figure;
plot(Wv, DiffuseToGlobalBOA);
title('Diffuse/Global Ratio, Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'DiffuseToGlobalIrradBOA', Rev, FileExts, TagFontProperties);

%% Read in the BWTek data for comparison 
% Uncertainty over the clock - is there a clock in the instrument ?
load(['..\Data\BWtekData\BWTekData' OverpassDate '.mat']);
figure;
plot(Wv, fastsmooth(GlobalBOAirrad, 30), BWTekDataOrdered(1).Wavelength,  BWTekDataOrdered(1).IrradiancemWcm2nm1/1000, ...
    BWTekDataOrdered(end-10).Wavelength,  BWTekDataOrdered(end-10).IrradiancemWcm2nm1/1000);
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
grid();
%SaveTaggedPlots(GitDescr, ResultsFolder,  'BWTekIrradBOA', Rev, FileExts, TagFontProperties);

%% Read in and plot the ASD data - definitely UTC
load(['..\Data\ASDIrrad\ASDIrradS3on' OverpassDate '.mat']);
figure;
plot(Wv, fastsmooth(GlobalBOAirrad, 120), ASDIrradMean.Wv, ASDIrradMean.RadData/10000); % Converting to W/cm^2/nm
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
xlim([350 1000])
legend('MODTRAN Smoothed', 'ASD', 'location', 'best');
grid();
%% Compute the water-leaving radiance
% First read in the R_rs values from Mark Matthews
RoodeRrsAll = dlmread(RrsFile, '\t');
RrsWv = RoodeRrsAll(:,1);
Rrs06 = RoodeRrsAll(:, RrsColumns);  % Select for day
% figure;
% plot(RrsWv, Rrs06);
% title(['R_{rs} for Roodeplaat Dam on ' OverpassDate]);
% xlabel('Wavelength [nm]');
% ylabel('R_{rs} [sr^{-1}]');
% legend('P1','P2','P3','P4', 'location', 'best')
% grid();

%% Interpolate R_rs to flux/sc7 wavelength grid and plot
Rrs06Interp = interp1(RrsWv, Rrs06, Wv, 'linear');
figure;
plot(Wv, Rrs06Interp);
title(['R_{rs} for Roodeplaat Dam on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('R_{rs} [sr^{-1}]');
legend('P1','P2','P3','P4', 'location', 'best')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'RrsASD', Rev, FileExts, TagFontProperties);

%% Compute/plot water-leaving radiance by multuplying Rrs by the GlobBOAirrad
Lw = repmat(GlobalBOAirrad, 1, size(Rrs06Interp, 2)) .* Rrs06Interp;
% plot water-leaving radiance
figure
plot(Wv, Lw);
title(['Water-leaving Radiance at BOA, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w [W/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'LwAtBOA', Rev, FileExts, TagFontProperties);

%% Compute/plot water-leaving radiance at TOA by multiplying by the path
% transmittance
LwTOA = repmat(S2.sc7.TRANS, 1, size(Lw, 2)) .* Lw;
figure;
plot(Wv, LwTOA);
title(['Water-leaving Radiance at TOA, Roodeplaat ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w at TOA [W/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'RemainingLwAtTOA', Rev, FileExts, TagFontProperties);

%% Compute the sky radiance as seen by reflection
% H1 = sensor, H2 = target
S2Sky = S2.AttachAlb(Roode1AA, 1); % NSURF = 1, area-averaged albedo
S2Sky.SetCaseName(['S2Rood' OverpassDate 'Sky']);
% Modify the path geometry to looking upwards using the ANGLE parameter
S2Sky.PHI = 0;  % Submit to ANGLE input
S2Sky.ANGLE = OZA;
SunRelAzToSky = OAA + 180 - SAA;
if SunRelAzToSky < 0
    SunRelAzToSky = SunRelAzToSky + 360;
end
if SunRelAzToSky > 360
    SunRelAzToSky = SunRelAzToSky - 360;
end
S2Sky.PARM1 = SunRelAzToSky; % Set the solar-relative viewing azimuth
S2Sky.Run;

% Interpolate the water reflectance to the wavelength grid
WaterReflRhoInterp = interp1(WaterReflRho(:,1), WaterReflRho(:,2), Wv, 'pchip');
WaterReflectedSkyRadiance = S2Sky.sc7.TOTALRAD .* WaterReflRhoInterp; % microwatts/sr/cm^2/nm

%% Plot the waterleaving radiance and water-reflected sky radiance together
% Remember to convert to common units of microwatts/sr/cm^2/nm
TotalRadianceBOA = Lw * 1e6 + repmat(WaterReflectedSkyRadiance, 1, size(Lw, 2));  % microwatts/sr/cm^2/nm
figure;
plot(Wv, Lw*1e6, Wv, WaterReflectedSkyRadiance);
title(['Water-leaving and Sky-Reflected Radiance at BOA, S2 on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('L_w and Sky-Reflected L at BOA [\muW/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'Sky-Reflected', 'location', 'best');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'LwAndWaterReflLAtBOA', Rev, FileExts, TagFontProperties);

%% Compute and plot total radiance at TOA
TotalLTOA = TotalRadianceBOA .* repmat(S2.sc7.TRANS, 1, size(Lw, 2)) ...
                    + repmat(S2.sc7.TOTALRAD, 1, size(Lw, 2));  % microwatts/sr/cm^2/nm
figure;                
plot(Wv, TotalLTOA);
title(['Total Radiance at TOA, S2 on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Total Radiance at TOA [\muW/sr/cm^2/nm]')
legend('P1','P2','P3','P4', 'location', 'best');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'TotalLatTOA', Rev, FileExts, TagFontProperties);


%% Plot LwTOA over TotalLTOA
figure;
plot(Wv, 1e6 * LwTOA ./ TotalLTOA)
title(['L_w Over Total L at TOA, S2 on ' OverpassDate])
xlabel('Wavelength [nm]');
ylabel('L_w at TOA / Total L')
legend('P1','P2','P3','P4', 'location', 'best');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'LwOverTotalLatTOA', Rev, FileExts, TagFontProperties);


%% Extract signals for S2 bands 1 to 10 (B1 ... B8, B8A, B9)


% Just do a loop rather that acrobatics with matrix dimensions
% Final Destination si ChanLTOA - channel radiance at TOA
ChanLTOA = zeros(size(S2FltInterp, 2), size(TotalLTOA, 2));
for iChan = 1:numel(S2FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = repmat(S2FltInterp(:,iChan), 1, size(TotalLTOA, 2));
    ChanLTOAProduct = ChanFilter .* TotalLTOA;
    ChanLTOA(iChan, :) = trapz(Wv, ChanLTOAProduct) / S2FltInterpIntegral(iChan);
end
% Convert to mW/sr/m^2/nm from microwatts/sr/cm^2/nm, a factor of 10
ChanLTOAmW = ChanLTOA * 10;

% Read some water dominated pixels from the S2 image
S2SNAPpixels = ReadSNAPpinData(WaterSNAPpixels, ...
    'all_refl', 'B([0-9]+[A]?)');
S2BandLegends = [S2SNAPpixels.all_refl_toks{1:10}];
S2SNAPpixels.all_refl = S2SNAPpixels.all_refl(:,1:10);  % Take only first 10 channels (VNIR)
ChanWv = [443 490 560 665 705 740 783 842 865 945];
%% Plot the S2 TOA radiances with MODTRAN TOA radiances
% First have to convert S2 reflectances to radiances.
figure;
plot(ChanWv, ChanLTOAmW, ChanWv, S2SNAPpixels.all_radiance', 'o');
title(['TOA Channel Radiance : S2 on ' OverpassDate ' at Roodeplaat']);
xlabel('Wavelength [nm]');
ylabel('Channel Radiance [mW/sr/m^2/nm]');
legend('MOD P1', 'MOD P2', 'MOD P3', 'MOD P4', 'S2');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'S2andMODTRANTotalLatTOA', Rev, FileExts, TagFontProperties);
%% Comparison of means
if size(S2SNAPpixels.all_radiance, 1) > 1  % Take mean over all pixels
    MeanS2Radiances = mean(S2SNAPpixels.all_radiance);
else
    MeanS2Radiances = S2SNAPpixels.all_radiance;
end
MeanChanLTOAmW = mean(ChanLTOAmW, 2);
%plot(ChanWv, MeanChanLTOAmW, 'o-', ChanWv, MeanS2Radiances, 'x-');

%% Plot percentage errors
figure;
plot(ChanWv, 100*2*(MeanS2Radiances-MeanChanLTOAmW')./(MeanS2Radiances+MeanChanLTOAmW'), 'o');
title(['Percentage Error : S2 vs MODTRAN at TOA on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Error (Difference over Mean) [%]');
grid;
SaveTaggedPlots(GitDescr, ResultsFolder,  'RelativeErrorMODTRANvsS2', Rev, FileExts, TagFontProperties)

%% Compare the solar flux data
GlobalTOAsolirrad = S2.flx.DirectSol(:,2);  % W/cm^2/nm
GlobalTOAsolirrad = repmat(GlobalTOAsolirrad, 1, size(S2FltInterp, 2)); % W/cm^2/nm
ChanGlobTOAsolirrad = trapz(Wv, GlobalTOAsolirrad .* S2FltInterp) ./ S2FltInterpIntegral;
ChanGlobTOAsolirrad = ChanGlobTOAsolirrad * 1000 * 10000;  % Convert to mW/m^2/nm as for S2 data
ChanGlobTOAsolDNI = ChanGlobTOAsolirrad ./ cos(deg2rad(SZA));   % Get the direct normal irradiance
% Now read in the S2 data
S2SolarIrrad = ReadSNAPpinData('..\Data\Sentinel2\S2SolarFluxDataAtRoodeplaatOn20160605.txt', 'all_solar_flux', 'solar_flux_band_([0-9]+)');
% Unfortunately not in order, so determine the order
BandOrder = str2double([S2SolarIrrad.all_solar_flux_toks{:}]);
S2SolarFlux = mean(S2SolarIrrad.all_solar_flux);
% And reorder
NewOrder = sortrows([BandOrder' [1:numel(BandOrder)]']);
S2SolarFlux = S2SolarFlux(NewOrder(:,2));
%% Plot ratio of MODTRAN solar DNI at TOA to S2 product.
figure
plot(ChanWv, ChanGlobTOAsolDNI./S2SolarFlux, 'o');
title(['Solar DNI at TOA : MODTRAN / S2 (Spectrum ' SolarSpectrum ' )']);
xlabel('Wavelength [nm]')
ylabel('Solar DNI MODTRAN / S2');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'TOASolarDNIFluxMODTRANoverS2', Rev, FileExts, TagFontProperties);

%% Comparison of means after correcting to S2 solar flux
% Percentage errors
MeanChanLTOAmWCorr = MeanChanLTOAmW' .* S2SolarFlux ./ ChanGlobTOAsolDNI;
figure;
plot(ChanWv, 100*2*(MeanS2Radiances-MeanChanLTOAmWCorr)./(MeanS2Radiances+MeanChanLTOAmWCorr), 'o');
title(['Flux-Corrected Percentage Error : S2 vs MODTRAN at TOA on ' OverpassDate]);
xlabel('Wavelength [nm]');
ylabel('Error (Difference over Mean) [%]');
grid;
SaveTaggedPlots(GitDescr, ResultsFolder,  'RelErrorMODTRANfluxcorrectedVsS2', Rev, FileExts, TagFontProperties);

%% Now attempt a retrieval of Lw at BOA
RadianceBOAWithoutLw = WaterReflectedSkyRadiance;  % microwatts/sr/cm^2/nm
LTOAWithoutLw = RadianceBOAWithoutLw .* S2.sc7.TRANS ...
                    + S2.sc7.TOTALRAD;  % microwatts/sr/cm^2/nm
%% Compute channel radiances for LTOA Without Lw
ChanLTOAWoLw = zeros(size(S2FltInterp, 2), size(LTOAWithoutLw, 2));
for iChan = 1:numel(S2FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = repmat(S2FltInterp(:,iChan), 1, size(LTOAWithoutLw, 2));
    ChanLTOAProductWoLw = ChanFilter .* LTOAWithoutLw;
    ChanLTOAWoLw(iChan, :) = trapz(Wv, ChanLTOAProductWoLw) / S2FltInterpIntegral(iChan);
end
% Convert to mW/sr/m^2/nm from microwatts/sr/cm^2/nm, a factor of 10
ChanLTOAmWWoLw = ChanLTOAWoLw * 10;

%% Compute and plot the percentage difference between S2 at TOA
% and MODTRAN at TOA without Lw
MeanChanLTOAmWWoLw = mean(ChanLTOAmWWoLw, 2);
figure
plot(ChanWv, MeanChanLTOAmWWoLw, 'o-', ChanWv, MeanS2Radiances, 'x-');
% Correct for solar flux differences
MeanChanLTOAmWWoLwCorr = MeanChanLTOAmWWoLw' .* S2SolarFlux ./ ChanGlobTOAsolDNI;
RetrievedLwAtTOA = MeanS2Radiances - MeanChanLTOAmWWoLwCorr;
figure;
plot(ChanWv, RetrievedLwAtTOA);

%% Now back to BOA by dividing by the path transmittance
% First need to average the path transmittance over each of the S2
% bands
S2BandPathTransmittance = zeros(size(S2FltInterp, 2), size(LTOAWithoutLw, 2));
for iChan = 1:numel(S2FltInterpIntegral)
    % Grab the channel flt function and replicate up to the number of
    % water sample positions
    ChanFilter = S2FltInterp(:,iChan);
    PathTransmittanceProduct = ChanFilter .* S2.sc7.TRANS;
    % Compute Weighted mean path transmittance
    S2BandPathTransmittance(iChan) = trapz(Wv, PathTransmittanceProduct) / S2FltInterpIntegral(iChan);
end
% Retrieve Lw at BOA
RetrievedLwAtBOA = RetrievedLwAtTOA ./ S2BandPathTransmittance'; % mW/sr/m^2/nm
%% Plot retrieved L_w at BOA
% Lw is in units of W/sr/cm^2/sr, so multiply by 1000 * 100 * 100
figure;
plot(Wv, Lw*100*100*1000, ChanWv, RetrievedLwAtBOA, 'ko-');  % in mW/sr/m^2/nm
title('S2 Retrieved L_w and Calculated L_w at BOA');
xlabel('Wavelength [nm]');
ylabel('L_w at BOA [mW/sr/m^2/nm]');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'S2RetrievedAndMeasuredLwAtBOA', Rev, FileExts, TagFontProperties);

print([ResultsFolder filesep 'S2RetrievedAndMeasuredLwAtBOARev' Rev '.pdf'], '-dpdf');
print([ResultsFolder filesep 'S2RetrievedAndMeasuredLwAtBOARev' Rev '.png'], '-dpng');


