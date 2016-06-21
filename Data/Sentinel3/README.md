# Sentinel 3 Data for Roodeplaat 2016 Campaign
Data related to Sentinel 3 is place in this folder.

Some important files here include the following :

## S3_0506_Roodeplaat_pinsMM.placemark
Pin placement for S3 image on 20160605 provided by Mark Matthews. Related email reads as follows:
Had a detailed look at the product. Looks like there is a row of duplicate pixels. Spatial offset is roughly 2 pixels to the right. Shapefile overlay clearly shows this. Verified on Hartbeespoort. 
Optimal pixel according to graphs is my pixel 1 or your pin 5. It has the lowest radiance, signifying it’s over the water. Nice colour scale clearly shows the water. Pin 1 and 4 are duplicates of Pin 2 and Pin 5. 
My suggestion: use only Pin 5 = Pin a. That is the lowest value and likely to be free (almost) of land contamination outside stray light.

The Pin 5 = Pin a pixel is the one used for the Rev001 computation of TOA radiance for the S3 overpass on 20160605.

## Sentinel3SRF2011Cam4.flt
MODTRAN "filter" file definitions for the Sentinel 3 spectral response functions.

