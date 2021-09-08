% This function detrends GPS carrier phase
% Inputs:
% phi --- Carrier phase time sequence
% fs --- Sampling frequency
% ds --- Standard deviation time interval
function [detPhi stdPhi] = butterdetrend(phi, fs, ds)
N = length(phi);
ButterFreq=.05; % Cutoff frequency in Hz.
ButterOrder=2; % Filter order.
[B,A]=butter(ButterOrder,ButterFreq/(fs/2),'high');
detPhi=filter(B,A,phi); detPhi=filter(B,A,detPhi); detPhi=filter(B,A,detPhi);
j = 0;
for i=1:(N-ds*fs)
j = j+1; temp = detPhi(i:i+ds*fs); stdPhi(j) = nanstd(temp);
end
end