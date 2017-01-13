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

%% Task 1: Determine Bandwidths of Interest

%%% Generate the Linear FM Signal (s(t)_Sim)
load lfm_sig.mat;

%%% Diagnostic Plot to Find Indices of Interest
% plot(10*log10(NormSpec(1:length(fplot))));
% set(gca,'YLim',[-20 2])
% grid
% xlabel ('Index')
% ylabel ('|IFFT|^{ 2} (dB)')

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

%% Task 2: Generate Realizations of AWGN

%%% Generate n_{AWGN}(t) and (S/N)_{Sim} with SNR in [-10.0 0.5 10.0] 
coeff_noise_db = -10:.5:10;
coeff_noise = 10.^(coeff_noise_db/10);

agwn = randn(1,length(LFMsig));
agwn = sqrt(coeff_noise*Pave)'*agwn;
sn_sim = agwn(1:end,:)+LFMsig;

%%% 2.A: Verify the SNR of (S/N)_{Sim}
power_ave_time_agwn = var(agwn,0,2) + mean(agwn,2).^2;

snr_sim = Pave./power_ave_time_agwn;
snr_sim_db = 10*log10(snr_sim);

%%% 2.B: Determine Bandwidth Containing ALL of the Simulated Noise Power
fft_pts = 2^nextpow2(length(agwn(1,:)));

psd_agwn = abs(fft(agwn,fft_pts,2)/length(agwn(1,:))).^2;

power_ave_freq_agwn = sum(psd_agwn,2)/2;

%% Task 3: Filter (S/N)_Sim

%%% 3.A: Generate a Butterworth LPF of Order 4 with BW matching W_{-3.0dB}

[butter_b butter_a] = butter(2,[f_3db_low,f_3db_high]/(FSamp/2));

% plot the filter response to verify functionality

[mag_butter, f_butter] = freqz(butter_b, butter_a, 5000, FSamp);

figure(1) % Question about how you want to show this
plot(f_butter,20*log10(abs(mag_butter))); %Question about 20 vs 10
set(gca,'YLim',[-5 2])
grid
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude Response of a 4^{th} order Butterworth BPF')

%%% 3.B: Filter the LFM Signal and (S/N)_{Sim} Realizations 

LFM_sig_filt = filtfilt(butter_b, butter_a, LFMsig);

for n=1:size(sn_sim,1)
        sn_sys(n,:) = filtfilt(butter_b, butter_a, sn_sim(n,:));
end

%%% 3.C: Plot (S/N)_{Sim} vs (S/N)_{Sys} in Time

figure(2)
plot(TimVec,sn_sim(21,:)); % 0.0 dB point
hold on
plot(TimVec,sn_sys(21,:));
grid
xlabel('Time (s)')
ylabel('Amplitude')
title('S/N_{Sim} vs. S/N_{Sys}')
legend(['S/N_{Sim}'; 'S/N_{Sys}'])
hold off

%%% 3.D: Plot (S/N)_{Sim} vs (S/N)_{Sys} in Freq

psd_sn_sim = abs(fft(sn_sim(21,:))/length(sn_sim(21,:))).^2;
psd_sn_sys = abs(fft(sn_sys(21,:))/length(sn_sys(21,:))).^2;

figure(3)
plot(fplot,10*log10(psd_sn_sim(1:length(fplot))));
hold on
plot(fplot,10*log10(psd_sn_sys(1:length(fplot))));
set(gca,'YLim',[-30 -20])
set(gca,'XLim',[1.2e6 5.2e6])
grid
xlabel ('Frequency (Hz)')
ylabel ('|IFFT|^{ 2} (dB)')
title('PSD of S/N_{Sim} vs. S/N_{Sys}')
legend(['S/N_{Sim}'; 'S/N_{Sys}'])
hold off

%%% 3.E: Plot LFM sig vs LFM sig filt

%% Task 4: Compare and Discuss the LFM time-domain at input and output

%%% 4.A Structure

%%% 4.B Power Characteristics

%%% 4.C Spectral Characteristics

%% Task 5: Compare and Discuss the time-domain noise at input and output

%%% 5.A Structure

%%% 5.B Power Characteristics

%%% 5.C Spectral Characteristics
