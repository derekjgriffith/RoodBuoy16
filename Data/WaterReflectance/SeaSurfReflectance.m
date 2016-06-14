function rho = SeaSurfReflectance(Phi_v, Theta_v, SZA, WindSpeed, InterpMethod)
% SeaSurfReflectance : Effective reflectance of the sea surface
%
% Usage :
%  >>> rho = SeaSurfReflectance(Phi_v, Theta_v, SZA, WindSpeed)
%      or
%  >>> rho = SeaSurfReflectance(Phi_v, Theta_v, SZA, WindSpeed, InterpMethod)
%
% Interpolates the effective reflectance of the sea surface using the
% data table published by Mobley :
%    Mobley, C.D., 2015.  Polarized Reflectance and Transmittance Properties of Wind-blown Sea Surfaces,
%    Applied Optics (in review)
%  NOTE: This is a preliminary table for use during review of the above paper.
%    These rho values are for fully developed seas, sun's direct rays parallel to the wind direction,
%    a Rayleigh, single-scattering sky at 550 nm, and wind speeds of 0, 2, 4, 5, 6, 8, 10, 12, 14, and 15 m/s.
%    The final table will be placed online once the above paper has been accepted for publication.
%  rho is computed from the I components of the Stokes radiance vectors: rho = I(surface reflected)/I(sky) [nondimensional]
%  Theta_v is the polar viewing angle measured from 0 at the nadir in degrees.
%  Phi_v is the azimuthal viewing angle measured from the sun at phi_v = 0
%  SZA is the solar zenith angle in degrees
%  WindSpeed is the wind speed in m/2
%
% The Mobley table can be downloaded from :
%   http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
%
% This function interpolates rho using interpn() on a grid of values where
% each of the inputs can be a row vector of values. 
%
% The optional input InterpMethod can be one of 'nearest', 'linear',
% 'spline' or 'cubic' as for interpn().
%
% Note that wavelength dependence is not implemented in this data/function.
% For more information on wavelength dependence, see the paper by Lee et
% al.
% Lee, Z. P., Y.-H. Ahn, C. D. Mobley, and R. Arnone, 2010. Removal of 
%   surface-reflected light for the measurement of remote-sensing 
%   reflectance from an above-surface platform. Optics Express, 
%   18(25), 26313-26324.

persistent SeaSurfRho;  % Load data only the first time the function is called.

if isempty(SeaSurfRho)
    load('MobleySeaSurfaceRho2015.mat');
end
if ~exist('InterpMethod', 'var')
    InterpMethod = 'linear';
end
if any(Phi_v > 180 | Phi_v < 0)
    warning('SeaSurfReflectance:PhiWarning', ...
        'Valid range for azimuthal angle Phi_v is 0 to 180 deg');
end
if any(Theta_v > 87.5 | Theta_v < 0)
    warning('SeaSurfReflectance:ThetaWarning', ...
        'Valid range for polar angle Theta_v is 0 to 87.5 deg');
end
if any(SZA > 87.5 | SZA < 0)
    warning('SeaSurfReflectance:SZAWarning', ...
        'Valid range for solar zenith angle SZA is 0 to 87.5 deg');
end
if any(WindSpeed > 15 | WindSpeed < 0)
    warning('SeaSurfReflectance:WindSpeedWarning', ...
        'Valid range for wind speed is 0 to 15 m/s');
end
rho = interpn(SeaSurfRho.Phi_v, SeaSurfRho.Theta_v, SeaSurfRho.SZA, SeaSurfRho.WindSpeed, SeaSurfRho.rho, ...
              Phi_v, shiftdim(Theta_v, 1), shiftdim(SZA, 2), shiftdim(WindSpeed, 3), InterpMethod);


end

