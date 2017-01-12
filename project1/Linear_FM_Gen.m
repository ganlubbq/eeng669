%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      LFM Signal Generator for EENG 699, Dig Comm I, Project #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Created: 22 Dec 04 (Dr. Michael A. Temple)
%  Last Checked/Updated: 23 Dec 16 (Dr. Michael A. Temple)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all

Tau = 90e-6;        % Signal Duration / Pulse Width (Sec)
CarFrq = 3250e3;     % Carrier/Center Frequency (Hz)
Pave = 1.58;        % Average LFM Signal Power (Watts)
FSamp = 60*CarFrq;  % Sample Frequency (Samp/Sec)
TSamp = 1/FSamp;       % Sample Time / Spacing (Sec)
TimVec = 0:TSamp:Tau-TSamp; % Time Vector

% Generate Simulated LFM Signal
pcr = 350;                           % Pulse Compression Ratio
flow = CarFrq-((1/Tau)/2)*pcr;       % Lower LFM Freq Parameter
fhgh = CarFrq+((1/Tau)/2)*pcr;        % Upper LFM Freq Parameter
coeff1 = (fhgh-flow)/(2*(Tau-TSamp));                            
coeff2 = flow;
LFMsig = sqrt(2*Pave)*cos((2*pi*(coeff1*TimVec.^2 + coeff2*TimVec)));

PowEstTim = var(LFMsig)+mean(LFMsig)^2; % Estimated Power (Time Domain)

% Generate and Plot LFM Signal: Time-Domain and PSD Responses

disp(' ')
disp('          !!!!! User Beware !!!!!')
disp(' Confirm/Verify Reliabiliity of Plotted Data')
disp(' ')

figure
subplot(3,1,1) % Time-Domain Waveform
plot(TimVec,LFMsig)
set(gca,'YLIM',[-3 3])
grid
title('LFM Signal Time Domain Response')
xlabel ('Time (Sec)')
ylabel ('Amp (v)')

LFMpsd = abs(fft(LFMsig)/length(LFMsig)).^2;   % PSD via |DFT|^2
MaxPSD = max(abs(LFMpsd));
MaxPSDdB = 10*log10(MaxPSD);
PowEstFrq = sum(LFMpsd);   % Estimated Power (Freq Domain)

subplot(3,1,2) % Un-Normalized Spectral PSD
delf = 1/Tau;    % Freq Sample Spacing
fplot=[0:delf:FSamp/2];    %
fplot(length(fplot))=[];
% Plot Norm PSD (Pos Freqs ONLY)
plot(fplot,10*log10(LFMpsd(1:length(fplot))));
set(gca,'YLim',[(MaxPSDdB-20) (MaxPSDdB+2)])
grid
title({['LFM Signal PSD: UN-NORMALIZED Continuous Plot'],...
    ['(Pos Freqs Only, Expanded About F_{Car} = ',num2str(CarFrq),' Hz)']})
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')

subplot(3,1,3) % Normalized Spectral PSD
delf = 1/Tau;    % Freq Sample Spacing
fplot=[0:delf:FSamp/2];    %
fplot(length(fplot))=[];
NormSpec = LFMpsd(1:length(fplot))/MaxPSD; % Normalize PSD
% Plot Norm PSD (Pos Freqs ONLY)
plot(fplot,10*log10(NormSpec(1:length(fplot))));
set(gca,'YLim',[-20 2])
grid
title({['LFM Signal PSD: NORMALIZED Continuous Plot '],...
    ['(Pos Freqs Only, Expanded About F_{Car} = ',num2str(CarFrq),' Hz)']})
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


