%% Read CEOS Thuiller Spectrum and resample for MODTRAN
% Must also convert irradiance from mW/m^2/nm to W/cm^2/cm^-1
clear all
close all
WvUnits = 'nm';
RadUnits = 'mW/m^2/nm';
Thui = dlmread('Solar_irradiance_Thuillier_2002.csv');
Wv = Thui(:,1);
ThuiRad = Thui(:,2);
plot(Wv, ThuiRad);
% Compute the wavenumber in cm^-1
WvNum = flipud(1e7./Wv);
% First convert from mW/m^2 to W/cm^2
ThuiRadWcm = flipud(ThuiRad/1000/100/100);
% Taken from ReadMODTRANfl7.m in MZDDEHg
% Notes: MODTRAN tape7 spectral radiance outputs are in units of W/cm^2/sr/cm^-1
%        i.e. given with respect to wavenumber. Conversion to W/cm^2/sr/nm (with respect to
%        wavelength) involves an important subtlety. Notice that the dimensionalty of the two units
%        is different, and this is the clue that a simple conversion factor does not apply.
%        if wn is wavenumber (in cm^-1) and wv is wavelenth (in nm), then wv = 1e7 wn^-1 and
%        delta_wv = -1e7 wn^-2 * delta_wn and delta_wn = -wn^2 / 1e7
%        Converting radiance to W/cm^2/sr/cm^-1 to W/cm^2/sr/nm therefore involves multiplying by
%        the factor wn^2/1e7 and flipping the result to get it into order of increasing wavelength.
%
%        Typical code is
%        >> [wavenumbers, wnradiance, h] = ReadMODTRANfl7('tape7', {'TOTAL RAD'});
%        >> wavelengths = flipud(1e7 ./ wavenumbers);
%        >> wvradiance = flipud(wnradiance .* wavenumbers.^2 ./ 1e7);
% MODTRAN file header must be as follows
%   FREQ     SOLAR IRRADIANCE
%  (CM-1)    (W CM-2 / CM-1)
ThuiRadWn = 1e7 * ThuiRadWcm ./ (WvNum.^2);
plot(WvNum, ThuiRadWn);
% Resample to 0.1 cm^-1 over the available range
WvNum_p1 = min(round(WvNum * 10)./10):0.1:max(WvNum);
%WvNum_p1 = min(WvNum):0.1:max(WvNum);
ThuiRadWn_p1 = interp1(WvNum, ThuiRadWn, WvNum_p1, 'linear', 'extrap');
plot(WvNum_p1, ThuiRadWn_p1);
fid = fopen('SUNp1CEOSThuillier2005.dat', 'wt');
fprintf(fid, '   FREQ     SOLAR IRRADIANCE\n');
fprintf(fid, '  (CM-1)    (W CM-2 / CM-1)\n');
fprintf(fid, '%8.2f     %11.4e\n', [WvNum_p1; ThuiRadWn_p1]);
fclose(fid);