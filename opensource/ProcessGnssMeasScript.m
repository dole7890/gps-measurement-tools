close all; clear all; clc
addpath('myfunction')
%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you: 
<<<<<<< Updated upstream
prFileName = 'pseudoranges_log_2016_06_30_21_26_07.txt'; %with duty cycling, no carrier phase
=======
% prFileName = 'gnss_log_2021_07_06_17_14_24.txt'; %with duty cycling, no carrier phase
prFileName = 'gnss_log_2021_07_06_17_13_55.txt'; %with duty cycling, no carrier phase
>>>>>>> Stashed changes
% prFileName = 'pseudoranges_log_2016_08_22_14_45_50.txt'; %no duty cycling, with carrier phase
% as follows
% 1) copy everything from GitHub google/gps-measurement-tools/ to 
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
<<<<<<< Updated upstream
dirName = '~/Documents/MATLAB/gpstools/opensource/demoFiles';
=======
dirName = 'C:\Users\Dong Kyeong Lee\Desktop\Research\ION2021\AndroidFaultMonitoring\210706_Baseline_Data\Antenna1\XiaomiBlack';
dirName = 'C:\Users\Dong Kyeong Lee\Desktop\Research\ION2021\AndroidFaultMonitoring\210706_Baseline_Data\Antenna3\XiaomiWhite';
>>>>>>> Stashed changes
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

%% data
%To add your own data:
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
param.llaTrueDegDegM = [37.55261226, 127.07380044, 87.577];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
if 0
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
save('temp.mat','gnssRaw','gnssAnalysis')
else
%     load('temp.mat')
    load('ant1_mi8.mat')
end
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
% [gnssMeas] = ProcessGnssMeas(gnssRaw);

for ii=1:0
%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph);

%% plot Pvt results
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
h5 = figure;
PlotPvtStates(gpsPvt,prFileName);
end

for ii=1:0
%% Plot Accumulated Delta Range 
if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
    [gnssMeas]= ProcessAdr(gnssMeas);
    h6 = figure;
<<<<<<< Updated upstream
    PlotAdr(gnssMeas,prFileName,colors);
    [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
=======
%     PlotAdr(gnssMeas,prFileName,colors);
    PlotAdr(gnssMeas,prFileName);
    [adrResid,sat_pos]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
>>>>>>> Stashed changes
    h7 = figure;
    PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
end
<<<<<<< Updated upstream
=======
end


for errormonitor = 1:1
[gnssMeas] = ProcessAdr(gnssMeas);
[adrResid,sat_pos]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
timeSeconds =adrResid.FctSeconds-adrResid.FctSeconds(1);%elapsed time in seconds

for ii = 1:length(adrResid.ResidM(1,:))/2
    if sum(~isnan(gnssMeas.AdrM(:,ii*2))) == 0
        continue
    end
    
    figure()
    % Detrended carrier
    plot(timeSeconds,adrResid.ResidM(:,ii*2-1))
    
    hold on
    
    % Android LOL Indicator
    iCs = find(bitand(gnssMeas.AdrState(:,ii*2-1),2^2));
%     plot(timeSeconds(iCs),adrResid.ResidM(iCs-1,ii*2-1).*ones(length(iCs),1),'xk')
    plot(timeSeconds(iCs),zeros(length(iCs),1),'xk')
    
    % Geometry Free
    GFC = diff(gnssMeas.AdrM(:,ii*2-1)-gnssMeas.AdrM(:,ii*2));
    idx3 = find(abs(GFC)>1);
    plot(timeSeconds(idx3),adrResid.ResidM(idx3,ii*2-1).*ones(length(idx3),1),'r*')
    
    % Butterworth
    idx=find(isnan(gnssMeas.AdrM(:,ii*2-1)));
    for jj=1:length(idx)
        if idx(jj)==1
            gnssMeas.AdrM(idx(jj),ii*2-1)=0;
            continue
        end
        gnssMeas.AdrM(idx(jj),ii*2-1)=gnssMeas.AdrM(idx(jj)-1,ii*2-1);
    end
    [detPhi,stdPhi] = butterdetrend(gnssMeas.AdrM(:,ii*2-1),1,60);
    idx2 = find(abs(detPhi)>0.3);
    idx0 = idx2(1);
    for kk=2:length(idx2)
        if idx2(kk)-idx2(kk-1) >100
            idx0 = [idx0;idx2(kk)];
        end
    end
%     plot(timeSeconds(idx0),adrResid.ResidM(idx0,ii*2-1).*ones(length(idx0),1),'go')
    plot(timeSeconds(idx0),zeros(length(idx0),1),'go')
    xlabel('Time (sec)')
    ylabel('Detrended Carrier (m)')
    if isempty(idx3)
        legend('Carrier','Flag: Android','Flag: Butterworth')
    else
        legend('Carrier','Flag: Android','Flag: GFC','Flag: Butterworth')
    end
    set(gca,'fontsize',12)
    grid on
    
%     figure();plot(detPhi)
%     continue
    % Code Rate Minus Carrier Rate
    figure()
    CMC = diff(gnssMeas.PrM(:,ii*2-1))-diff(gnssMeas.AdrM(:,ii*2-1));
    idx = find(isfinite(CMC));
    plot(abs(CMC),'b');hold on
    threshold = (4*gnssMeas.PrSigmaM(:,ii*2-1).^2+4*gnssMeas.AdrSigmaM(:,ii*2-1).^2).^0.5;
    plot(timeSeconds,threshold,'r')
    ylim([0 20])
    
    xlabel('Time (sec)')
    ylabel('Code Rate - Carrier Rate(m/s)')
    grid on
    set(gca,'fontsize',12)
    
    % Doppler Minus Carrier Rate
    figure()
    DMC = gnssMeas.PrrMps(2:end,ii*2-1)-diff(gnssMeas.AdrM(:,ii*2-1));
%     idx = find(isfinite(DMC));
    DMC = DMC - DMC(idx(1));
    plot(DMC,'b');hold on
    threshold = (2*gnssMeas.PrrSigmaMps(:,ii*2-1).^2+4*gnssMeas.AdrSigmaM(:,ii*2-1).^2).^0.5;
    plot(timeSeconds,threshold,'r')
    ylim([0 1e-2])
    
    xlabel('Time (sec)')
    ylabel('Doppler - Carrier (m/s)')
    grid on
    set(gca,'fontsize',12)
    
    
end
end

save('ant3_mi8.mat','gnssRaw','gnssAnalysis','gnssMeas','sat_pos')


>>>>>>> Stashed changes
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
