%% Project #1: Simulated Filter Characterization
%
% This M-file is written to characterize the performace of basic MATLAB
% filtering operations and assessing their impact on simulated signals
% corrupted with AWGN.
%
% Inputs:
%   N/A
%
% Outputs:
%
% Dependancies:
%   Linear_FM_Gen.m
%
% Author:  2d Lt Taylor Bodin
% Class: Digital Communications I, EENG 669, Winter 2017
close all; clear; clc;
rng(20170111); %seed the RNG for reproducable results
set(0,'defaultfigureWindowStyle','docked')

%% Task 1: Determine Bandwidths of Interest

%%% Generate the Linear FM Signal (s(t)_Sim)
load lfm_sig.mat;

%%% 1.A.i: -3.0 dB Bandwidth Based on Center Frequency Power

% Setup Variables
index_CarFrq = floor(CarFrq/delf) + 1; % Carrier frequency is the middle of the signal
power_CarFrq = 10*log10(NormSpec(index_CarFrq));

% find the 3.0 and 10.0 dB down levels
down_3db = power_CarFrq - 3.0; % -4.3526 dB
down_10db = power_CarFrq - 10.0; % -11.3526 dB

% Graphically determined interpolation points for each level
lfm_low_3db = 10*log10(NormSpec(123:124)); 
lfm_high_3db = 10*log10(NormSpec(463:464));

lfm_low_10db = 10*log10(NormSpec(112:113)); 
lfm_high_10db = 10*log10(NormSpec(475:476));

% 3 dB
f_3db_low = interp1(lfm_low_3db, fplot(123:124),down_3db);
f_3db_high = interp1(lfm_high_3db, fplot(463:464),down_3db);
bw_3db = f_3db_high - f_3db_low;

% 10 dB
f_10db_low = interp1(lfm_low_10db, fplot(112:113),down_10db);
f_10db_high = interp1(lfm_high_10db, fplot(475:476),down_10db);
bw_10db = f_10db_high - f_10db_low;

%%% 1.C: Bandwidth Containing ALL of the Simulated LFM signal power
bw_all = fplot(end)-fplot(1);

%%% 1.D: Average LFM Signal Power Contained in W_{-3.0dB} and W_{-10.0dB}

% W_3db Average Power
index_3db_low = 123;
index_3db_high = 464;
power_3db_bw = sum(LFMpsd(index_3db_low:index_3db_high));

% W_10db Average Power
index_10db_low = 112;
index_10db_high = 476;
power_10db_bw = sum(LFMpsd(index_10db_low:index_10db_high));

%% Task 2: Generate (S/N)_{sim}

%%% 2.A: LFMsig

% Time Domain Imported from lfm_sig.mat
% - LMFsig

% Spectral Domain Imported from lfm_sig.mat
% - LFMpsd

% Power measures imported from lfm_sig.mat
% - PowEstFrq
% - PowEstTim
% - Pave

%%% 2.B: Generate AGWN

% Coefficients for SNR in [-10.0 0.5 10.0] 
coeff_noise_db = -10:.5:10;
coeff_noise = 10.^(coeff_noise_db/10);

% Time Domain
agwn = randn(1,length(LFMsig));
agwn = sqrt(coeff_noise*Pave)'*agwn;

% Spectral Domain
fft_pts = 2^nextpow2(length(agwn(1,:)))/2; % Used for all subsequent ffts
psd_agwn = abs(fft(agwn,fft_pts,2)/length(agwn(1,:))).^2;

% Power
power_ave_time_agwn = var(agwn,0,2) + mean(agwn,2).^2;
power_ave_freq_agwn = sum(psd_agwn,2)/2;

%%% 2.C: Generate LFMsig_sim

% Time Domain
LFMsig_sim = agwn(1:end,:)+LFMsig;

% Spectral Domain
for n=1:size(LFMsig_sim,1)
        psd_LFMsig_sim(n,:) = ...
            abs(fft(LFMsig_sim(n,:))/size(LFMsig_sim,2)).^2;
end

% Power
power_ave_time_LFMsig_sim = var(LFMsig_sim,0,2) + mean(LFMsig_sim,2).^2;
power_ave_freq_LFMsig_sim = sum(psd_LFMsig_sim,2);

%%% 2.D: Calculate SNR_{sim}
snr_sim = Pave./power_ave_time_agwn;
snr_sim_db = 10*log10(snr_sim);

%% Task 3: Generate (S/N)_{sys}

%%% 3.A: Generate a Butterworth LPF with BW matching W_{-3.0dB}
[butter_b butter_a] = butter(4,[f_3db_low,f_3db_high]/(FSamp/2));

% plot the filter response to verify functionality
% note: filt filt produces a slightly different effective filter
[mag_butter, f_butter] = freqz(butter_b, butter_a, 5000, FSamp);

figure(1) % Question about how you want to show this
plot(f_butter,10*log10(abs(mag_butter).^2)); %Question about 20 vs 10
set(gca,'YLim',[-5 2])
line([f_3db_low f_3db_low],get(gca,'YLim'),'Color',[1 0 0])
line([f_3db_high f_3db_high],get(gca,'YLim'),'Color',[1 0 0])
line(get(gca,'Xlim'),[-3 -3],'Color',[1 0 0])
text(f_3db_low, -3.2, [num2str(f_3db_low,'%10.3e') ' Hz'], 'FontSize', 12);
text(f_3db_high, -3.2, [num2str(f_3db_high,'%10.3e') ' Hz'], 'FontSize', 12);
grid
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude Response of the butter(4,W_n) Butterworth BPF')

%%% 3.B: Filter LFMsig

% Time Domain
LFMsig_filt = filtfilt(butter_b, butter_a, LFMsig);

% Spectral Domain
psd_LFMsig_filt = abs(fft(LFMsig_filt,fft_pts,2)/length(LFMsig_filt(1,:))).^2;

% Power

%%% 3.C: Filter AWGN

% Time Domain
for n=1:size(agwn,1)
        agwn_sys(n,:) = filtfilt(butter_b, butter_a, agwn(n,:));
end

% Spectral Domain
psd_agwn_filt = abs(fft(agwn_sys,fft_pts,2)/length(agwn_sys(1,:))).^2;

% Power
power_ave_time_sys_agwn = var(agwn_sys,0,2) + mean(agwn_sys,2).^2;

%%% 3.D: Filter the AWGN + LFMsig

% Time Domain done in loop in 3.B
for n=1:size(LFMsig_sim,1)
        LFMsig_sys(n,:) = filtfilt(butter_b, butter_a, LFMsig_sim(n,:));
end

% Spectral
psd_sn_sys = abs(fft(LFMsig_sys(21,:))/length(LFMsig_sys(21,:))).^2;

% Power
power_ave_time_sys_lfm = var(LFMsig_filt,0,2) + mean(LFMsig_filt,2).^2;

%%% 3.E: Calculate SNR_{sys}
snr_sys = power_ave_time_sys_lfm./power_ave_time_sys_agwn;
snr_sys_db = 10*log10(snr_sys);


%% Task 4: Compare and Discuss the time-domain LFMsig at input and output

%%% 4.A Structure
figure(5)

h(1) = subplot(1,2,1);
plot(TimVec,LFMsig);
grid
xlabel ('Time (s)')
ylabel ('Amplitude')
title('LFM Signal')

h(2) = subplot(1,2,2);
plot(TimVec,LFMsig_filt)
grid
xlabel ('Time (s)')
ylabel ('Amplitude')
title('Filtered LFM Signal')

%%% 4.D Spectral Characteristics
figure(6)

h(1) = subplot(1,2,1);
plot(fplot,10*log10(LFMpsd(1:length(fplot))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('PSD of the LFM Signal')

h(2) = subplot(1,2,2);
plot(fplot,10*log10(psd_LFMsig_filt(1:length(fplot)))); %/max(abs(psd_LFMsig_filt)
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('PSD of the Filtered LFM Signal')

%linkaxes(h,'x')
xlim(h,[1e6 5.5e6])
ylim(h,[-40 -25])

%%% 4.C Power Characteristics

%% Task 5: Compare and Discuss the time-domain noise at input and output

%%% 5.A Structure
figure(7)

h(1) = subplot(1,2,1);
plot(TimVec,agwn(21,:));
grid
xlabel ('Time (s)')
ylabel ('Amplitude')
title('Noise Structure Prior to Filtering')

h(2) = subplot(1,2,2);
plot(TimVec,agwn_sys(21,:))
grid
xlabel ('Time (s)')
ylabel ('Amplitude')
title('Noise Structure After Filtering')

%%% 5.B Spectral Characteristics
figure(8)
h(1) = subplot(1,2,1);
plot(fplot,10*log10(psd_agwn(21,1:length(fplot))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('Noise PSD Prior to Filtering')

h(2) = subplot(1,2,2);
plot(fplot,10*log10(psd_agwn_filt(21,1:length(fplot))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('Noise PSD After Filtering')

linkaxes(h,'x')
xlim(h,[0 1e7])

%%% 5.C Power Characteristics

%% Task 6: Compare and Discuss LFMsig + AWGN

%%% 6.A: Structure
figure(2)

h(1) = subplot(1,3,1);
plot(TimVec,LFMsig); % 0.0 dB point
grid
xlabel('Time (s)')
ylabel('Amplitude')
title('The LFM Signal')

h(2) = subplot(1,3,2);
plot(TimVec,LFMsig_sim(21,:)); % 0.0 dB point
grid
xlabel('Time (s)')
ylabel('Amplitude')
title('The LFM Signal in the Simulated Channel')

h(3) = subplot(1,3,3);
plot(TimVec,LFMsig_sys(21,:));
grid
xlabel('Time (s)')
ylabel('Amplitude')
title('The LFM Signal Filtered in the Simulated System')

linkaxes(h,'y')
ylim(h(1),[-6 6])

%%% 6.B: Spectral
figure(3)

h(1) = subplot(1,3,1);
plot(fplot,10*log10(NormSpec(1:length(fplot))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('Normalized PSD of the LFM Signal')

h(2) = subplot(1,3,2);
plot(fplot,10*log10(psd_LFMsig_sim(21,1:length(fplot)) ...
    ./max(abs(psd_LFMsig_sim(21,:)))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('Normalized PSD of the LFM Signal in the Simulated Channel')

h(3) = subplot(1,3,3);
plot(fplot,10*log10(psd_sn_sys(1:length(fplot))./max(abs(psd_sn_sys))));
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('Normalized PSD of the LFM Signal Filtered in the Simulated System')

linkaxes(h,'y')
xlim(h,[0 7e6])
ylim(h,[-15 0])

%%% 6.C: Power
figure(4)

plot(snr_sim_db,snr_sys_db);
set(gca,'XLim',[-10 10])
set(gca,'YLim',[3 25])
grid
xlabel ('Simulated SNR (dB)')
ylabel ('System SNR (dB)')
title('SNR_{Sim} vs. SNR_{Sys}')