%% BS antenna pattern
% Pablo Caballero Garcés
% 30/03/15

function [ gain_sectorized ] = antenna_gain(angle)
gain_sectorized(:)=-min(12.*(([60 180 300]-angle)./70).^2,20);

%s1=60;
%gain_sectorized(1)=-min(12*((s1-angle)/70)^2,20);

%s2=180;
%gain_sectorized(2)=-min(12*((s2-angle)/70)^2,20);

%s3=300;
%gain_sectorized(3)=-min(12*((s3-angle)/70)^2,20);

end
