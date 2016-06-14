%% Read the sea surface reflectance factors by Mobley
% The header in the file reads as follows :
%  Table of rho values created with FFT surfaces and polarized ray tracing as described in
%    Mobley, C.D., 2015.  Polarized Reflectance and Transmittance Properties of Wind-blown Sea Surfaces,
%    Applied Optics (in review)
%  NOTE: This is a preliminary table for use during review of the above paper.
%    These rho values are for fully developed seas, sun's direct rays parallel to the wind direction,
%    a Rayleigh, single-scattering sky at 550 nm, and wind speeds of 0, 2, 4, 5, 6, 8, 10, 12, 14, and 15 m/s.
%    The final table will be placed online once the above paper has been accepted for publication.
%  rho is computed from the I components of the Stokes radiance vectors: rho = I(surface reflected)/I(sky) [nondimensional]
%  Theta_v is the polar viewing angle measured from 0 at the nadir in degrees.
%  Phi_v is the azimuthal viewing angle measured from the sun at phi_v = 0
%  in degrees.
%  Looping order is
%     wind speed
%        sun zenith angle
%           theta_v
%              phi_v
%  Data blocks of 118 theta_v, phi_v, rho(theta_v, phi_v) records on an (f7.1,f9.1,e15.4) format are separated by 
%    records giving the wind speed and sun zenith angle on an (12x,f5.1,20x,f5.1) format.
% -------------------------------
%  Theta_v    Phi_v       rho
% -------------------------------
% WIND SPEED =  0.0; SUN ZENITH ANGLE =  0.0
fid = fopen('WaterRhoTable_AO2015.txt', 'rt');
WindSpeed = [];
SZA = [];
allrho = [];
while ~feof(fid)
    line = fgetl(fid);
    tok = strtok(line);
    if strcmp(tok, 'WIND') % Found block of data
        % obtain wind sped and solar zenith angle
        WsSZA = sscanf(line, 'WIND SPEED =  %f; SUN ZENITH ANGLE = %f');
        WindSpeed = [WindSpeed WsSZA(1)];
        SZA = [SZA WsSZA(2)];
        Data = textscan(fid, '%f %f %f');
        Theta_v = Data{1};
        Phi_v = Data{2};
        rho = Data{3};
        % Have to append copies of the first value of rho
        rho = [repmat(rho(1), numel(unique(Phi_v)) - 1, 1); rho];
        allrho = [allrho; rho];
    end
end
fclose(fid);

WindSpeed = unique(WindSpeed);
SZA = unique(SZA);
Theta_v = unique(Theta_v);
Phi_v = unique(Phi_v);

% allrho1 = reshape(allrho, [numel(WindSpeed), numel(SZA), numel(Theta_v), numel(Phi_v)]);
% The innermost loop variable occurs first in indexing 
allrho2 = reshape(allrho, [numel(Phi_v), numel(Theta_v), numel(SZA), numel(WindSpeed)]);
rho = allrho2;
% Save the data for use in a function
SeaSurfRho.rho = rho;
SeaSurfRho.Theta_v = Theta_v;
SeaSurfRho.Phi_v = Phi_v;
SeaSurfRho.SZA = SZA;
SeaSurfRho.WindSpeed = WindSpeed;
SeaSurfRho.Ref = 'Mobley, C.D., 2015.  Polarized Reflectance and Transmittance Properties of Wind-blown Sea Surfaces, Applied Optics';
SeaSurfRho.AxisDefinitions = {'Theta_v is the polar viewing angle measured from 0 at the nadir in degrees', ...
    'Phi_v is the azimuthal viewing angle measured from the sun at phi_v = 0','Solar Zenith Angle','Wind Speed'};
SeaSurfRho.AxisOrder = {'Phi_v','Theta_v','SZA','WindSpeed'};
SeaSurfRho.AxisUnits = {'deg','deg','deg','m/s'};
save MobleySeaSurfaceRho2015.mat SeaSurfRho

