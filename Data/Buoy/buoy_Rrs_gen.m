%Buoy_Rrs_Gen
%Draft script to pull in Luke generated .csv buoy data for Roodeplaat
%Deployment Jun2016

%Define wavelength
lambda=[320:5:895];

%Get file name
	[csv_name, in_pname]=uigetfile('*.csv', 'Choose a buoy csv file');

%Load rads data
All = csvread(csv_name,1,0)

%Assign into acquisitions
Rad0945(:,1)=All(:,2)
Rad0945(:,2)=All(:,3)
Rad0945(:,3)=All(:,4)

Rad1000(:,1)=All(:,5)
Rad1000(:,2)=All(:,6)
Rad1000(:,3)=All(:,7)

Rad1015(:,1)=All(:,8)
Rad1015(:,2)=All(:,9)
Rad1015(:,3)=All(:,10)

%Convert from Trios units 
Rad0945=Rad0945./10;
Rad1000=Rad1000./10;
Rad1015=Rad1015./10;

%Ed(find(T2_data(:,116)>1),:)=[];
%Lu_1(find(T2_data(:,116)>1),:)=[];
%Lu_2(find(T2_data(:,116)>1),:)=[];
%T_meta(find(T2_data(:,116)>1),:)=[];

%Dirty R calc
%R_ss=T2_data./T1_data;
%Remove strange spectra with high NIR
%R_ss(find(R_ss(:,116)>0.005),:)=[];

%Set depth factors
z1=0.4;
delta_z=0.5;


%Calculate Ku data asssuming delta_Z=0.4
Ku0945(:,1)=(log(Rad0945(:,2))-log(Rad0945(:,3)))/(delta_z);
Ku1000(:,1)=(log(Rad1000(:,2))-log(Rad1000(:,3)))/(delta_z);
Ku1015(:,1)=(log(Rad1015(:,2))-log(Rad1015(:,3)))/(delta_z);


%Calculate Rrs data
Rrs(:,1)=(Rad0945(:,2).*exp((Ku0945.*(0.4)).*(0.98./1.334.^2)))./Rad0945(:,1);
Rrs(:,2)=(Rad1000(:,2).*exp((Ku1000.*(0.4)).*(0.98./1.334.^2)))./Rad1000(:,1);
Rrs(:,3)=(Rad1015(:,2).*exp((Ku1015.*(0.4)).*(0.98./1.334.^2)))./Rad1015(:,1);

%Calculate Chl from Algal0/709 algorithm
C709=237.5.*exp(-2.13.*(Rrs(70,:)./Rrs(79,:)));

%Figure of all spectral data
figure
   
subplot(2,2,1);
set(gca,'FontName','times','FontSize',12);
plot(lambda,Rad0945(:,1),'b'); hold; grid
plot(lambda,Rad1000(:,1),'g');
plot(lambda,Rad1015(:,1),'r');
axis([300 900 0 100])
ylabel('\muW cm^{-2} nm^{-1}')
xlabel('\lambda [nm]')
title('Irradiance Ed');

subplot(2,2,2)
set(gca,'FontName','times','FontSize',12);
plot(lambda,Rad0945(:,2:3),'b'); hold; grid
plot(lambda,Rad1000(:,2:3),'g')
plot(lambda,Rad1015(:,2:3),'r')
axis([300 900 0 1])
ylabel('\muW cm^{-2} sr^{-1} nm^{-1}')
xlabel('\lambda [nm]')
%legend('Lu z1','Lu z2');
title('Radiance Lu');

subplot(2,2,3)
set(gca,'FontName','times','FontSize',12);
plot(lambda,Ku0945(:,1),'b'); hold; grid
plot(lambda,Ku1000(:,1),'g');
plot(lambda,Ku1015(:,1),'r');
axis([300 900 0 5])
ylabel('m^{-1}')
xlabel('\lambda [nm]')
title('K_{Lu}');

subplot(2,2,4)
set(gca,'FontName','times','FontSize',12);
plot(lambda,Rrs(:,1),'b'); hold; grid
plot(lambda,Rrs(:,2),'g');
plot(lambda,Rrs(:,3),'r');
axis([300 900 0 0.015])
ylabel('sr^{-1}')
xlabel('\lambda [nm]')
title('Reflectance R_{rs}');

  