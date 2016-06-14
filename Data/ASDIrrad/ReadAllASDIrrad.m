Directory = 'D:\Projects\Calval\Roodeplaat\Data\ASDIrrad\20160605\';
ASDIrrad = ASD([Directory 'Irrad*.asd.irr.txt']);
% Print the dates and times of measurements (Computer clock should be set
% to UTC). The overpass time is 07:42:31
%% Find measurements near the overpass time
disp(strvcat(ASDIrrad.DateTime));
ASDSerialDates = [ASDIrrad.SerialDate];
% Want to find measurements within 3 minutes of the overpass
OverPassSerialDate = datenum([2016 6 5 7 42 33]);
MinSerialDate = datenum([2016 6 5 7 39 33]);
MaxSerialDate = datenum([2016 6 5 7 45 33]);
iMeasure = find(ASDSerialDates > MinSerialDate & ASDSerialDates < MaxSerialDate);
ASDIrradNearOverpass = ASDIrrad(iMeasure);
ASDIrradNearOverpass.Plot();
ASDIrradMean = ASDIrradNearOverpass.Mean;
ASDIrradMean.Plot();
save ASDIrradS3on20160605.mat ASDIrradMean