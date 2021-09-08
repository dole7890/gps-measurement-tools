clc; clear all; close all;

load('ant1_mi8.mat')
% [gnssMeas] = ProcessGnssMeas(gnssRaw);
A = gnssMeas;
sat_pos_A = sat_pos;
load('ant3_mi8.mat')
% [gnssMeas] = ProcessGnssMeas(gnssRaw);
B = gnssMeas;
sat_pos_B = sat_pos;

t1 = A.FctSeconds;
t2 = B.FctSeconds;
t = intersect(t1,t2);

for ii=1:length(t)
    idx1(ii) = find(t1==t(ii));
    idx2(ii) = find(t2==t(ii));
end

sv1 = 2;
sv2 = 12;

ADR_Ap = A.AdrM(idx1,sv1*2-1);
ADR_Bp = B.AdrM(idx2,sv1*2-1);

ADR_Aq = A.AdrM(idx1,sv2*2-1);
ADR_Bq = B.AdrM(idx2,sv2*2-1);

sat_pos_A = sat_pos_A(idx1,:);

DD = ADR_Ap - ADR_Bp - ADR_Aq + ADR_Bq;

res=[];
for ii=1:length(t)
x1 = sat_pos_A(ii,(sv1*6-4):(sv1*6-2));
x2 = sat_pos_A(ii,(sv2*6-4):(sv2*6-2));
xa = [-3052145.347399311, 4039497.553048760, 3866248.475357282];
xb = [-3052145.008637840, 4039497.592240163, 3866248.700323119];

G_1 = [(x1(1)-xa(1));...
       (x1(2)-xa(2));...
       (x1(3)-xa(3))];
G_2 = [(x2(1)-xa(1));...
       (x2(2)-xa(2));...
       (x2(3)-xa(3))];

geo_range_p_a(ii) = norm(G_2);

G_1 = G_1/norm(G_1);
G_2 = G_2/norm(G_2);
G = G_2 - G_1;

res(ii) = DD(ii) - G'*(xb-xa)';
end
idx = find(isfinite(res));
% res = res - nanmean(res);%res(idx(1));
res = res - res(idx(1));
figure()
plot(res)

return
figure()
for ii=1:length(gnssMeas.Svid)/2
    idx = find(isfinite(gnssMeas.AdrM(:,(ii-1)*2+1)));
    plot(gnssMeas.FctSeconds(idx),ones(length(idx),1)*gnssMeas.Svid((ii-1)*2+1))
    hold on
end
return
for ii=1:length(A.PrM(1,:))/2
    
    
    figure()
    
    
    return
    A.AdrM(:,ii*2-1)
    B.AdrM(:,ii*2-1)
end
% 
% function svXyzTrxM = satpos(gnssMeas,
%     svid = gnssMeas.Svid;
% 
%     ttxSeconds = gnssMeas.tTxSeconds(i,j);
%     if isnan(ttxSeconds)
%         continue %skip to next
%     end
%     [gpsEph,iSv]= ClosestGpsEph(allGpsEph,svid(j),gnssMeas.FctSeconds(i));
%     if isempty(iSv)
%         continue; %skip to next
%     end
%     %compute pr for this sv
%     dtsv = GpsEph2Dtsv(gpsEph,ttxSeconds);
%     ttxSeconds = ttxSeconds - dtsv;%subtract dtsv from sv time to get true gps time
% 
%     %calculate satellite position at ttx:
%     [svXyzTtxM,dtsv]=GpsEph2Xyz(gpsEph,[weekNum(i),ttxSeconds]);
%     %in ECEF coordinates at trx:
%     dtflightSeconds = norm(xyz0M - svXyzTtxM)/GpsConstants.LIGHTSPEED;
%     svXyzTrxM = FlightTimeCorrection(svXyzTtxM, dtflightSeconds);
% end