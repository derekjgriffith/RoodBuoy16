% Run the script for the particular measurement and conditions before running
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
MODStartWv = 350; % Start wavelength in MODTRAN
MODStopWv = 1000; % Stop wavelength in MODTRAN, includes up to S3 b20
MODDV = 0.05; % Spectral interval in MODTRAN
MODFWHM = 0.1; % Run at 0.1 nm smoothing (convolution) in MODTRAN
MODWv = MODStartWv:MODDV:MODStopWv; % Anticipated MODTRAN output wavelengths
Mod5.ParallelFriendly(true);
% This case is visible/near-infared (VIS/NIR) wavelengths
Diff2Glob = Mod5;    % Get a completely empty case instance
% Set up name and short description
Diff2Glob = Diff2Glob.SetCaseName(['Diff2GlobRood' MeasureDate]); % The SetCaseName method is the only way to set the CaseName property
Diff2Glob.CaseDescr = ['Diff2Glob Measure at Roodeplaat Dam on ' MeasureDate];

% Plotting sizes
theFontSize = 14;
theLineWidth = 2.0;

% Note that if a card is required, ALL parameters on that card must be set,
% even if the parameters are not used.

% Set up Card 1 (mandatory - main radiative transport)
Diff2Glob.MODTRN = 'M';     % MODTRAN band model
Diff2Glob.SPEED = 'S';      % Slow algorithm
Diff2Glob.BINARY = 'f';     % Output will be ASCII
Diff2Glob.LYMOLC = ' ';     % Exclude 16 auxiliary trace gases
Diff2Glob.MODEL = 3;        % Mid-latitude winter canned atmosphere 
Diff2Glob.ITYPE = 3;        % Slant path to ground
Diff2Glob.IEMSCT = 4;       % Solar radiation only, no thermal scatter
Diff2Glob.IMULT = -1;       % Include multiple scatter, computed at H2 (target/pixel)
Diff2Glob.M1 = 0;           % Temperature/pressure default to MODEL (Mid-latitude winter profile)
Diff2Glob.M2 = 0;           % Water vapor defaults to MODEL profile
Diff2Glob.M3 = 0;           % Ozone defaults to MODEL profile
Diff2Glob.M4 = 0;           % Methane defaults to MODEL profile
Diff2Glob.M5 = 0;           % Nitrous oxide defaults to MODEL profile
Diff2Glob.M6 = 0;           % Carbon monoxide defaults to MODEL profile
Diff2Glob.MDEF = 0;         % Default O2, NO, SO2, NO2, NH3, and HNO3 species profiles.
Diff2Glob.I_RD2C = 0;       % Normal program operation - no user input for profiles
Diff2Glob.NOPRNT = 0;       % Minimize printing to Tape6 output file
Diff2Glob.TPTEMP = 0;       % Temperature at H2 - not important, only VIS/NIR
Diff2Glob.SURREF = '0.0';   % Set surface reflectance across whole spectrum
% Note that the setting of Card 1 parameters can be accomplished in a single
% call to the Set method as follows :
% Diff2Glob = Diff2Glob.Set('MODTRN', 'M','SPEED', 'S', 'BINARY', 'f', 'LYMOLC', ' ','MODEL', 3, 'ITYPE', 3, 'IEMSCT', 2, ...
%                 'IMULT', 0, 'M1', 0, 'M2', 0, 'M3', 0, 'M4', 0, 'M5', 0 , 'M6', 0, ...
%                 'MDEF', 0, 'I_RD2C', 0, 'NOPRNT', 1, 'TPTEMP', 0, 'SURREF', '0.5');


% Set up Card 1A (mandatory - main radiative transport continued)
Diff2Glob.DIS = 't';        % Using DISORT multiple scattering algorithm
Diff2Glob.DISAZM = 't';     % Therefore also not using azimuth dependence in DISORT
Diff2Glob.DISALB = 'f';     % Don't calculate atmospheric correction data
Diff2Glob.NSTR = NSTR;        % Number of streams for DISORT
Diff2Glob.SFWHM = 0;        % Default solar irradiance data
Diff2Glob.CO2MX = 380;      % CO2 mixing ratio, xxx ppm by volume
Diff2Glob.H2OSTR = H2OSTR;  % Scale/set of water vapor amount
Diff2Glob.O3STR = O3STR;    % Scale/set of ozone profile amount
Diff2Glob.C_PROF = '0';     % No scaling of default molecular species profiles
Diff2Glob.LSUNFL = LSUNFL;     % Read specified solar irradiance data
Diff2Glob.LBMNAM = 'f';     % Don't read alternative band model file
Diff2Glob.LFLTNM = 't';     % Must read filter file specified
Diff2Glob.H2OAER = 'f';     % Don't bother to modify aerosol properties on the basis of H2OSTR
Diff2Glob.CDTDIR = 'f';     % Data files are in the default location
Diff2Glob.SOLCON = -1;      % Unity scaling of TOA solar irradiance, but apply seasonal correction
Diff2Glob.CDASTM = CDASTM;     % Angstrom law perturbation flag.
% Set up Angstrom law perturbations if provided
if exist('ASTMC', 'var')
    Diff2Glob.ASTMC = ASTMC; 
end
if exist('ASTMX', 'var')
    Diff2Glob.ASTMX = ASTMX;
end
if exist('ASTMO', 'var')
    Diff2Glob.ASTMO = ASTMO;
end
% Diff2Glob.AERRH
Diff2Glob.NSSALB = NSSALB;       % Manipulate single scattering albedo if required

% Deal with card 1A1 - user-defined solar spectrum
if strcmpi(LSUNFL, 'T')
    Diff2Glob.USRSUN = ['..' filesep 'DATA' filesep USRSUN];   
end

% Card 1B for single scattering albedo
if NSSALB > 0
    Diff2Glob.AWAVLN = AWAVLN;   % Wavelengths for single scattering albedo in microns
    Diff2Glob.ASSALB = ASSALB;   % Single scattering albdeo
end

% Set up Card 2 (mandatory - main aerosol and cloud options)
Diff2Glob.APLUS = '  ';     % Don't use flexible aerosol manipulations
Diff2Glob.IHAZE = IHAZE;        % Rural aerosol model, visibility = 23 km (modified below)
Diff2Glob.CNOVAM = ' ';     % Don't invoke NOVAM
Diff2Glob.ISEASN = 0;       % Use default seasonal aerosol tweaking
Diff2Glob.ARUSS = '   ';    % Don't use extended user-defined aerosol facility
Diff2Glob.IVULCN = 0;       % Background stratospheric aerosol profile
Diff2Glob.ICSTL = 1;        % Continental influence of maritime aerosols - not applicable to this case
Diff2Glob.ICLD = 0;         % No clouds or rain
Diff2Glob.IVSA = 0;         % Don't use Army Vertical Structure Algorithm for boundary layer aerosols
Diff2Glob.VIS = -AOT550;    % Negative of the AOT at 550
Diff2Glob.WSS = 2;          % Use default wind speed for named MODEL
Diff2Glob.WHH = 2;          % Use default 24 hr average wind speed for named MODEL
Diff2Glob.RAINRT = 0;       % Rain rate is zero (mm/hour), anyway no cloud/rain (ICLD)
Diff2Glob.GNDALT = GNDALT;       % Target surface (H2) is at sea level

% Set up Card 3 (mandatory - Line of sight geometry)
% To define path (LOS) geometry in this case use PHI, H1 and H2 (combination 3c in manual)
% Here we want a vertical path from ground to space - set everything to
% zero
Diff2Glob.H1 = 0;           % Not used in this case - we are using a slant path to ground/space
Diff2Glob.H2 = 0;           % km. Target pixel is at sea level
Diff2Glob.ANGLE = 0;        % Not used in this case. (Zenith angle at H1)
Diff2Glob.RANGE = 0;        % Not used in this case. Path length.
Diff2Glob.BETA = 0;         % Not used in this case. Earth centre angle.
Diff2Glob.RO = 0;           % Not used in this case. Radius of the Earth, will default to a reasonable value.
Diff2Glob.LENN = 0;         % Not used in this case. Short path/long path switch.
Diff2Glob.PHI = 0;        % degrees. Zenith angle at H2 (pixel/target) to H1 (satellite camera)

% Set up Card 3A1 (Solar scattering geometry, required for IEMSCT = 2)
Diff2Glob.IPARM = 12;       % Will specify relative solar azimuth angle and solar zenith angle below (PARM1 and PARM2)
Diff2Glob.IPH = 2;          % Use Mie-generated internal database for aerosol phase functions (???????????)
Diff2Glob.IDAY = MeasureDateVec(1:3); % Compute day number corresponding to day of Measure to adjust TOA solar irradiance
Diff2Glob.ISOURC = 0;       % The Sun is the extraterrestrial source of scattered radiation

% Set up Card 3A2 (Solar scattering geometry, also required for IEMSCT = 2)
SunRelAz = SAA - OAA;
if SunRelAz < 0
    SunRelAz = SunRelAz + 360;
end
Diff2Glob.PARM1 = SunRelAz;       % deg. The Sun azimuth is x deg (north through east positive) of LOS azimuth (H2 to H1)
Diff2Glob.PARM2 = SZA;      % deg. Sun zenith angle at H2 (target/pixel).
Diff2Glob.PARM3 = 0;        % Not used in this case.
Diff2Glob.PARM4 = 0;        % Not used in this case.
Diff2Glob.TIME = 0;         % Not used in this case.
Diff2Glob.PSIPO = 0;        % Not used in this case.
Diff2Glob.ANGLEM = 0;       % Not used in this case.
Diff2Glob.G = 0;            % Not used in this case. (Henyey-Greenstein asymmetry parameter)
% An alternative way of setting up cards 3A1 and 3A2 (shortwave source scattering geometry)
% is to use the method SetScatGeom defined as follows:
%   MC = MC.SetScatGeom(IPARM, IDAY, ISOURC, PARM, TIME, PSIPO, ANGLEM)
% In this case, the call would be :
%   Diff2Glob = Diff2Glob.SetScatGeom(12, [2009 11 2], 0, [50 20]);
% (TIME, PSIPO and ANGLEM are not used and will be set 0 by SetScatGeom)
% Note that IPH must still be set explicitly, as well as G.

% Set up Card 4 (mandatory - spectral range and resolution)
Diff2Glob.V1 = MODStartWv;         % Start of spectral computation range in nm (see FLAGS(1))
Diff2Glob.V2 = MODStopWv;         % End of spectral computation range in nm
Diff2Glob.DV = MODDV;         % Spectral increment in nm
Diff2Glob.FWHM = MODFWHM;         % Convolution filter width in nm
Diff2Glob.YFLAG = ' ';      % Not going to generate .plt or .psc files
Diff2Glob.XFLAG = ' ';      % Not going to generate .plt or .psc files
Diff2Glob.FLAGS(1) = 'N';   % Use nanometres for spectral units (FLAGS(1)).
Diff2Glob.FLAGS(4) = 'A';   % Put ALL radiance components in convolved data (tp7)
% Want to output fluxes as well
Diff2Glob.FLAGS(7) = 't';  % Output fluxes
Diff2Glob.MLFLX = 1;  % Only output fluxes at BOA and TOA
% An alternative way of setting up (most of) Card 4 is to use the SetSpectralRange method
% The call looks as follows:
%   MC = MC.SetSpectralRange(V1, V2, DV, FWHM, Units, ConvShape, FWHMisRelative)
% and in this case would look as follows :
%   Diff2Glob = Diff2Glob.SetSpectralRange(350, 650, 0.1, 2, 'N');
% (ConvShape and FWHMisRelative are not required)
% Note that it will still be necessary to set FLAGS(4) as this is not done
% by SetSpectralRange. YFLAG, XFLAG and DLIMIT can be set using the
% SetPlot method.

% Read area-averaged albedo retrieved from S3 on 2016-06-05 at Roodeplaat
S3Rood220160605RetrievedAlbedo = Mod5.ReadAlb('..\Data\ASDIrrad\SunPopRefl0605\S3Rood220160605.alb');
% Attach the albedo data, specifically the area-averaged albedo
Diff2Glob = Diff2Glob.AttachAlb(S3Rood220160605RetrievedAlbedo, 1); 

% Set up Card 5 (mandatory - Repeat option)
Diff2Glob.IRPT = 0;         % Stop program, only one sub-case in this run

% Now run the case (execute MODTRAN on the case)
Diff2Glob = Diff2Glob.Run;
% Record the optical depth
Diff2GlobTransTotalOD = Diff2Glob.sc7.DEPTH;
% Create a copy
Diff2GlobTrans = Diff2Glob;
% Switch off aerosols and run again
Diff2GlobTrans.IHAZE = 0;
Diff2GlobTrans = Diff2GlobTrans.Run;
Diff2GlobTransNoAerosolOD = Diff2GlobTrans.sc7.DEPTH;
Diff2GlobAOD = Diff2GlobTransTotalOD - Diff2GlobTransNoAerosolOD;
AODWv = Diff2GlobTrans.sc7.WAVLNM;

%% Plot the AOT vs measurements
figure;
plot(AODWv, Diff2GlobAOD, AOTwv, AOT, 'or', 'LineWidth', theLineWidth);
set(gca, 'LineWidth', theLineWidth);
set(gca, 'FontSize', theFontSize);
title(['Vertical AOD, Roodeplaat ' MeasureDate])
xlabel('Wavelength [nm]');
ylabel('AOD')
legend(['MODTRAN ' IHAZEModel], 'MicroTOPS', 'location', 'best');
axis([350, 900, 0, 1])
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'AOD', Rev, FileExts, TagFontProperties);
% Save AOT Measurements
save(['RoodeplaatAOT' MeasureDate IHAZEModel '.mat'], 'AODWv', 'Diff2GlobAOD', 'AOTwv', 'AOT');

%% Calculate and plot the downwelling spectral irradiance
Wv = Diff2Glob.flx.Spectral;  % nm
% Compute total downwelling
GlobalBOAirrad = Diff2Glob.flx.DownDiff(:,1) + Diff2Glob.flx.DirectSol(:,1);
figure;
plot(Wv, GlobalBOAirrad);
title('Total Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Total Irradiance [W/cm^2/nm]');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'TotalDownIrradBOA', Rev, FileExts, TagFontProperties);
% !!! Note that irradiances are in W/cm^2/nm, while
% Diff2Glob.sc7 radiances are in microwatts/sr/cm^2/nm

%% Also plot diffuse to global
DiffuseToGlobalBOA = Diff2Glob.flx.DownDiff(:,1) ./ GlobalBOAirrad;
figure;
plot(Wv, DiffuseToGlobalBOA);
title('Diffuse/Global Ratio, Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio');
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'DiffuseToGlobalIrradBOA', Rev, FileExts, TagFontProperties);

%% Read in ASD and BWTek diffuse to global and plot together
% ASD data saved by Saved by ReadPlotSunPopRefl.m
% save IrradDiffuseGlobalRatio20160605.mat IrradRatio IrradRatioMean
load(ASDSunPopDiff2GlobFile);

% BWTek data saved by ReadAllBWTekIrradExp0605.m
% save BWTekDataIrradExp20160605.mat BWTekDataUnordered Wavelength GlobMeanIrrad DiffMeanIrrad Diff2GlobRatio Diff2GlobRatioSmooth
load(BWTekSunPopDiff2GlobFile);

% Compute the misfit at the AOT control wavelengths
% Interpolate the ASD and MODTRAN results at the AOT wavelengths
Diff2GlobMOD = interp1(Wv, DiffuseToGlobalBOA, AOTwv, 'linear');
Diff2GlobASD = interp1(IrradRatioMean.Wv, IrradRatioMean.RadData, AOTwv, 'linear');
Diff2GlobBWTek = interp1(Wavelength, Diff2GlobRatioSmooth, AOTwv, 'linear');

Misfit = sqrt(sum((Diff2GlobMOD - Diff2GlobASD).^2))
BWTekMisfit = sqrt(sum((Diff2GlobMOD - Diff2GlobBWTek).^2))

figure;
plot(Wv, DiffuseToGlobalBOA, IrradRatioMean.Wv, IrradRatioMean.RadData, Wavelength, Diff2GlobRatioSmooth, ...
    AOTwv, Diff2GlobMOD, 'bo', AOTwv, Diff2GlobASD, 'go');
title('Diffuse/Global Ratio, Downwelling Irradiance at BOA');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio');
axis([350 1000 0 0.8]);
legend('MODTRAN', ['ASD Misfit ' num2str(Misfit)], 'BWTek', 'location', 'northeast')
grid();
SaveTaggedPlots(GitDescr, ResultsFolder,  'Diff2GlobMODvsASDvsBWTek', Rev, FileExts, TagFontProperties);

% Compute the misfit at the AOT control wavelengths
% Interpolate the ASD and MODTRAN results at the AOT wavelengths
Diff2GlobMOD = interp1(Wv, DiffuseToGlobalBOA, AOTwv, 'linear');
Diff2GlobASD = interp1(IrradRatioMean.Wv, IrradRatioMean.RadData, AOTwv, 'linear');

%% Plot the difference over the mean
% First have to sample to same wqavelength grid
Diff2GlobASDatMODWv = interp1(IrradRatioMean.Wv, IrradRatioMean.RadData, Wv);
figure;
Diffuse2GlobMODSmooth = fastsmooth(DiffuseToGlobalBOA, 200, 1, 1);
plot(Wv, 200 * (Diff2GlobASDatMODWv - Diffuse2GlobMODSmooth) ...
           ./ (Diff2GlobASDatMODWv + Diffuse2GlobMODSmooth));
%% Plot again with smoothed MODTRAN data       

figure;
plot(Wv, Diffuse2GlobMODSmooth, '--', IrradRatioMean.Wv, IrradRatioMean.RadData, Wavelength, Diff2GlobRatioSmooth, ':', ...
    AOTwv, Diff2GlobMOD, 'bo', AOTwv, Diff2GlobASD, 'go', AOTwv, Diff2GlobBWTek, 'rx', 'LineWidth', theLineWidth);
set(gca, 'LineWidth', theLineWidth);
set(gca, 'FontSize', theFontSize);
title('Diffuse/Global Ratio, Downwelling Irradiance');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio');
axis([350 1000 0.1 0.8]);
legend(['MODTRAN ' IHAZEModel], ['ASD, Misfit ' num2str(Misfit)], ['BWTek, Misfit ' num2str(BWTekMisfit)], 'location', 'northeast')
grid();
%SaveTaggedPlots('', ResultsFolder,  'DiffuseToGlobalIrradBOASmooth', Rev, FileExts, TagFontProperties);

%% Also do a percentage error plot

Diff2GlobASDatModWv = interp1(IrradRatioMean.Wv, IrradRatioMean.RadData, Wv, 'linear');
Diff2GlobBWTekAtModWv = interp1(Wavelength, Diff2GlobRatioSmooth, Wv, 'linear');
Diff2GlobErrorMODvsASD = 200 * (Diff2GlobASDatModWv - Diffuse2GlobMODSmooth) ...
           ./ (Diff2GlobASDatModWv + Diffuse2GlobMODSmooth);
Diff2GlobErrorMODvsBWTek = 200 * (Diff2GlobBWTekAtModWv - Diffuse2GlobMODSmooth) ...
           ./ (Diff2GlobBWTekAtModWv + Diffuse2GlobMODSmooth);
figure;       
plot(Wv, Diff2GlobErrorMODvsASD, 'b--', Wv, Diff2GlobErrorMODvsBWTek, 'g', 'LineWidth', theLineWidth);       
set(gca, 'LineWidth', theLineWidth);
set(gca, 'FontSize', theFontSize);
title('Diffuse/Global Ratio Error');
xlabel('Wavelength [nm]');
ylabel('Diffuse/Global Ratio Error [%]');
legend(['MODTRAN ' IHAZEModel ' vs ASD'], ['MODTRAN ' IHAZEModel ' vs BWTek']);
xlim([350 900]);
grid();

