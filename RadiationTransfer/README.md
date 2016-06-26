# Atmospheric Radiative Transfer
This folder contains scripts and data for atmospheric radiative transfer calculations related to the Roodeplaat/Sentinel campaign 2016.

Both MODTRAN and libRadtran have been used for this work. In order to use the MODTRAN scripts here, a copy of MODTRAN 5 will be required. See http://www.modtran5.com

For Matlab scripts which use MODTRAN, the Matlab/MODTRAN wrapper is also required, which can be found at https://github.com/derekjgriffith/matlab-modtran-5

The general procedure for taking in-situ measurements to TOA using MODTRAN is a s follows:

MODTRAN was used to propagate matchup Rrs measurements at Roodeplaat to TOA for specific overpasses of S2 and S3 (2016-06-05 and 2016-06-06).
Aerosol loading was moderate on these particular days (~0.5 at 500 nm), dominated presumably by biomass-burning with urban-transport-industrial component.
Aerosol Optical Thickness (AOT) and water-vapour column were measured at Roodeplaat with MicroTOPS sunphotometers, cross-calibrated with Aeronet Cimel CE318 at CSIR (~ 15 km away). 
Using AOD and water vapour measured in-situ, the area-averaged (~ 2km diameter) surface reflectance was retrieved from the S2/S3 image using MODTRAN. Found to have rather apparent errors due to mismatch in MODTRAN solar TOA spectrum with reference to solar flux provided in S3 product (even for Thuillier-derived spectrum in MODTRAN).
Area-averaged surface reflectance was used to calculate total atmospheric radiance at TOA for completely black pixel (black water, no water-surface reflectance). Also used this MODTRAN run to compute total downwelling irradiance at BOA.
Used total downwelling irradiance at BOA to compute water-leaving radiance Lw at BOA via Rrs measurements performed in-situ around the time of overpasses (MM to elucidate protocol).
Added water-surface reflected sky radiance (using an upward-looking MODTRAN run in the correct viewing geometry) to Lw to obtain total upwelling radiance above water at BOA. Water-surface reflectance was interpolated from the Mobley 2015 tables (assumed spectrally invariant, a point which may require addressing).
Multiplied total above-water upwelling radiance at BOA by atmospheric path transmittance to obtain water-target radiance at TOA. Added total atmospheric radiance to obtain total water-target apparent radiance (sensor-observed) at TOA.
All computations up to this point performed at full MODTRAN 5 spectral resolution.
Applied S2/S3 spectral response functions to compute band radiances and compared to S2 (via reflectance product) or S3 (directly) measurements at carefully selected pixels in the product.
In the case of S3, corrections were applied to align the MODTRAN solar flux to the solar flux reported in the S3 product.

The retrieval of Lw at BOA from the S2/S3 images was performed using exactly the same process as above except that the Lw component was left out of the computation (but not the water surface-reflected component). The band radiance shortfall at TOA between MODTRAN and S2/S3 measurement was than taken to be the Lw component at TOA. This was divided by the atmospheric path (band averaged) transmittance to obtain S2/S3 band spectral water-leaving radiance at BOA and compared to Lw at BOA computed via Rrs measurements.

For the near-nadir (~14 deg Observation Zenith Angle) overpass of S3 on 5 June, the canned MODTRAN urban aerosol model anchored by measured AOT at 550 nm was found to be quite competent at retrieving Lw at BOA in visible spectrum. However, the S3 overpass on 6 June had OZA of 50.6 deg and retrieval errors were much larger. Conditions also not as stable with a bit of cloud appearing between S3 and S2 overpass. It was necessary to enhance the MODTRAN aerosol model with Angstrom parameter and spectral Single Scattering Albedo (SSA) data to obtain a reasonable retrieval. SSA data was taken from Aeronet at CSIR (~ 15 km away), assumed to be representative of fairly large area. Further enhancements of the aerosol model may yield still better retrieval of Lw at BOA.

libRadtran can be obtained from http://www.libradtran.org

For Python scripts or Jupyter notebooks which use libRadtran, the Python MORTICIA package is required, which can be obtained at https://github.com/derekjgriffith/MORTICIA
