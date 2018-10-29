function [ tbs ] = cqi2tbs(cqi)
%3GPP TSG RAN1 #52bis,R1-081638,Shenzhen, China
%Annex 1

%% CQI to MCS

%mcss=[0 1 2 4 6 8 11 13 15 18 20 22 24 26 28]

%mcs=mcss(cqi+1)
%% MCS to TBS
%http://www.etsi.org/deliver/etsi_ts/136200_136299/136213/08.08.00_60/ts_136213v080800p.pdf
% 50 RBS
%tbss=[0 1 2 4 6 8 10 12 14 17 19 20 22 24 26]
%tbs=tbss(cqi+1)
Nprb=[0 1384 1384 2216 3624 ...
      5160 6968 8760 11448 ...
      14112 18336 21384 22920 ...
      27376 30576 36696];
tbs=Nprb(cqi+1)*1000; %this is bps

