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

%%% 1.A.i: -3.0 dB Bandwidth Based on Center Frequency Power
index_CarFrq = floor(CarFrq/delf) + 1; % Carrier frequency is the middle of the signal
power_CarFrq = 10*log10(NormSpec(index_CarFrq));

down_3db_cent = power_CarFrq - 3.0; % -4.3526 dB
f_3db_low_cent = 1.359e6; % graphically determined
f_3db_high_cent = 5.144e6;
bw_3db_cent = f_3db_high_cent - f_3db_low_cent;

%%% 1.A.ii: -3.0 dB Bandwidth Based on Max Power
down_3db_max = -3; % Max Power is 0.0 dB in NormSpec by definition
f_3db_low_max = 1.384e6;
f_3db_high_max = 5.151e6;
bw_3db_max = f_3db_high_max - f_3db_low_max;

%%% 1.B.i: -10.0 dB Bandwidth Based on Center Frequency Power
down_10db_cent = power_CarFrq - 10.0; % -11.3526 dB
f_10db_low_cent = 1.244e6; 
f_10db_high_cent = 5.267e6; 
bw_10db_cent = f_10db_high_cent - f_10db_low_cent;

%%% 1.B.ii: -10.0 dB Bandwidth based on Max Power
down_10db_max = -10;
f_10db_low_max = 1.267e6; 
f_10db_high_max = 5.244e6; 
bw_10db_max = f_10db_high_max - f_10db_low_max;

%%% 1.C: Bandwidth Containing ALL of the Simulated LFM signal power
bw_all = fplot(end)-fplot(1);

%%% 1.D: Average LFM Signal Power Contained in W_{-3.0dB} and W_{-10.0dB}

% Diagnostic Plot to Find Indices for Summing Up Average Power
% plot(10*log10(NormSpec(1:length(fplot))));
% set(gca,'YLim',[-20 2])
% grid
% xlabel ('Index')
% ylabel ('|IFFT|^{ 2} (dB)')

% W_3db Average Power
index_3db_low_max = 125;
index_3db_high_max = 462;
power_3db_bw = sum(LFMpsd(index_3db_low_max:index_3db_high_max));

% W_10db Average Power
index_10db_low_max = 114;
index_10db_high_max = 473;
power_10db_bw = sum(LFMpsd(index_10db_low_max:index_10db_high_max));

%% Task 2: Generate Realizations of AWGN

%%% Generate n_{AWGN}(t) and (S/N)_{Sim} with SNR in [-10.0 0.5 10.0] 
coeff_noise_db = -10:.5:10;
coeff_noise = 10^(coeff_noise_db/10);

awgn = randn(1,length(LFMsig));
awgn = sqrt(coeff_noise)'*awgn;

sn_sim = awgn(1:end,:)+LFMsig;

%%% 2.A: Verify the SNR of (S/N)_{Sim}
power_ave_time_awgn = var(awgn,0,2) + mean(awgn,2).^2;
power_ave_freq_awgn = var(awgn,0,2) + mean(awgn,2).^2;

snr_sim = power_ave_time_awgn./Pave;

%%% 2.B: Determine Bandwidth Containing ALL of the Simulated Noise Power

% Generate PSDs for awgn and sum over all of it see that it matches time and
% theoretical

%% Task 3: Filter (S/N)_Sim

%%% 3.A: Generate a Butterworth LPF of Order 4 with BW matching W_{-3.0dB}

bpf_butter = butter(4,[1.384e6,5.151e6]/(samplefreq/2));
% plot the filter response to verify functionality

%%% 3.B: Filter the LFM Signal and (S/N)_{Sim} Realizations 

%%% 3.C: Plot (S/N)_{Sim} vs (S/N)_{Sys}

%% Task 4: Compare and Discuss the LFM time-domain at input and output

%%% 4.A Structure

%%% 4.B Power Characteristics

%%% 4.C Spectral Characteristics

%% Task 5: Compare and Discuss the time-domain noise at input and output

%%% 5.A Structure

%%% 5.B Power Characteristics

%%% 5.C Spectral Characteristics
