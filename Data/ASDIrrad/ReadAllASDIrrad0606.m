Directory = '.\20160606\';
ASDIrrad = ASD([Directory 'Irrad*.asd.irr.txt']);
% Print the dates and times of measurements (Computer clock should be set
% to UTC). The overpass time is 07:16:36
%% Find measurements near the overpass time
disp(strvcat(ASDIrrad.DateTime));
ASDSerialDates = [ASDIrrad.SerialDate];
% Want to find measurements within 3 minutes of the overpass
OverPassSerialDate = datenum([2016 06 06 7 16 36]);
MinSerialDate = datenum([2016 6 6 7 14 36]);
MaxSerialDate = datenum([2016 6 6 7 18 36]);
iMeasure = find(ASDSerialDates > MinSerialDate & ASDSerialDates < MaxSerialDate);
ASDIrradNearOverpass = ASDIrrad(iMeasure);
ASDIrradNearOverpass.Plot();
ASDIrradMean = ASDIrradNearOverpass.Mean;
ASDIrradMean.Plot();
save ASDIrradS3on20160606.mat ASDIrradMean