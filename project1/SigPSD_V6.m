% =========================================================================
% =========================== SigPSD_V6 Function ==========================
% =========================================================================
%
%   function [SigPsd,FScale,SigPow] = SigPSD_V6(SigIn,FSamp,PSDFlr,...
%                                               SAve,NAve,SNorm,SPlot) 
%
% =========================================================================
%
%   Calculates Power Spectral Density (PSD) as the Magnitude of Fast
%   Fourier Transform (FFT) Coefficients Squared.  Resultant PSD (SigPsd),
%   Frequency Scale (FScale) and Average Power Est (SigPow) are output.
%
%   Created:  Apr 06, Dr. Michael A. Temple (MAT)
%   Modified: 19 Apr 11, MAT ... Enable Matrix Input/Output
%   Modified: 7 Sep 11, MAT ... Minor Change to Plotting Output
%   Modified: 30 Apr 13, MAT ... V3 ('PSDFlr' Input Added)
%   Modified: 21 Nov 14, MAT ... V4
%               Input/Ouput Compatibility with V3 Maintained
%               Updated SigIn Vector/Matrix Processing
%               Verified Freq Scales for Real/Complex (See Note Below) 
%               Output Plots Cleaned Up
%   Modified: 24 Aug 15, MAT
%               Complex Plots Cleaned Up
%   Modified: 1 Mar 16, MAT [V5]
%               Updated for Matlab 2015b Compatibility
%               E.g., FigDex Initialization
%   Modified: 1 Aug 16, MAT [V6]
%               Added 'conj' For Complex SigIn
%               Output for Complex SigIn Validated with TstSigs
%
%   NOTE:  Outputs Verified/Consistent with Matlab's 'spectrum' Generation 
%          Process for BOTH Complex-Valued and Real-Valued Input Signals
%          (21 Nov 14).  See/Run "Complex SigIn" Example #3 Below.
%
% =========================================================================
%                                   Inputs
% =========================================================================
%   SigIn   NSig x NSamp Matrix:
%               NSig = # Input Signals
%               NSamp = # Time Samples per Signal
%   FSamp   Sample Frequency (Samples/Sec)
%   PSDFlr  Output PSD Floor: Min[SigPsd] Value(s) = PSDFlr (dB)
%            0 --> No Flooring Applied [V5]
%           XX --> PSD Floor Set to -XX dB Below Max PSD Value.
%                  +/- XX Accepted as Input [V5]
%   SAve    PSD Averaging Control Variable (Sliding Window Average)
%           (1 = Average Applied, 0 = No Average Applied)
%   NAve    Number of Points Averaged Across
%           (Odd Number Required if Save = 1, Else, Doesn't Matter)
%   SNorm   Normalization Control Variable
%           (1 = Normalized, 0 = UnNormalized)
%   SPlot   Plot Control Variable (1 = Plot, 0 = No Plot).  Three plots
%           generated for SPlot=1:  A) Time Domain Response (Amplitude for
%           SigIn Real and Magnitude for SigIn Complex), B) SigIn PSD
%           (Spanning -FSamp/2 < f < FSamp/2 for Real-Valued SigIn and
%           0 < f < FSamp for SigIn Complex), and C) Angle of SigIn Fourier
%           Transform
%
% =========================================================================
%                                  Outputs
% =========================================================================
%   SigPsd  PSD of Input Signal Calc as |Fourier{SigIn}|^2
%   FScale  Resultant Frequency Scale for Output PSD
%   SigPow  Estimated Power Estimate Calc as Sum of PSD Components
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Example 1:  Complex SigIn with a D.C. Term and Sinusoids @ FCar1, FCar2
%               and FCar3 Hz.  No Averaging, UnNormalized, with Plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FSamp=25e4;
% XDelta=1/FSamp;
% NSamp=49500;
% TimVec=[0:XDelta:(NSamp-1)*XDelta];
% FCar1=33500;
% FCar2=41200;
% FCar3=63000;
% % FCar4 Intentionally Set Above Nyquist:  1) Appears INCORRECTLY as
% % Folder-Over Response at (FSamp/2)-12500 = 11250 Hz in Real-Valued PSD,
% % and 2) Appears CORRECTLY at f = 137500 in Complex PSD.
% FCar4=(FSamp/2)+12500; % Greater Than Nyquist
% SPow=2;
% %
% Sig1=(2*sqrt(SPow))*exp(-1i*2*pi*FCar1*TimVec);
% Sig2=(2*sqrt(SPow/4))*exp(-1i*2*pi*FCar2*TimVec);
% Sig3=(2*sqrt(SPow))*exp(-1i*2*pi*FCar3*TimVec);
% Sig4=(2*sqrt(SPow))*exp(-1i*2*pi*FCar4*TimVec);
% TSig=(SPow/4)+Sig1+Sig2+Sig3+Sig4;
%
% [TSigPSD,FScale,TSigPow]=SigPSD_V6(TSig,FSamp,40,0,0,0,1);
% [TSigPSD,FScale,TSigPow]=SigPSD_V6(real(TSig),FSamp,40,0,0,0,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Example 2:  Complex Linear Frequency Modulated (LFM) Signal.
%               No Averaging, Normalized, with Plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tau = 70e-3;    % Signal Duration / Pulse Width
% PCR = 210;      % Pulse Compression Ratio
% Pave = 2.37;    % Average LFM Signal Power
% FCntr=17500      % Approx Cntr Frq
% FLow = FCntr-((1/Tau)/2)*PCR;  % Lower LFM Freq
% FHgh = FCntr+((1/Tau)/2)*PCR;   % Upper LFM Freq
% Coeff1 = (FHgh-FLow)/(2*(Tau-XDelta));
% Coeff2 = FLow;
% Phz = 2*pi*(Coeff1*TimVec.^2 + Coeff2*TimVec);
% TSig = sqrt(2*Pave)*exp(-1i*Phz);
% [LFMPSD1,FScale1,LFMPow1]=SigPSD_V6(TSig,FSamp,40,0,0,0,1);
% [LFMPSD2,FScale2,LFMPow2]=SigPSD_V6(real(TSig),FSamp,40,0,0,0,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Example 3:  Complex SigIn with Multiple Sinusoids @ 1000, 7400, 18000,
%               43000 and (73000 > FSamp/2) Hz. No Averaging, UnNormalized,
%               with Plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       FSamp = 584000; % (2*73000)*4 = 4xNyquist
%       XDelta=1/FSamp;
%       NSamp=50000;
%       t = [0:XDelta:XDelta*(NSamp-1)];
%       TSig = 1.0*exp(1i*2*pi*1000*t) + 2.0* exp(1i*2*pi*7400*t) ...
%           + 1.0* exp(1i*2*pi*18000*t) + 2.0* exp(1i*2*pi*43000*t)...
%           + 2.0* exp(1i*2*pi*73000*t);
%       % Complex-Valued Sig PSD ... No Freq Fold-Over of f = 73000
%       [CSigPSD,CFScale,CSigPow]=SigPSD_V4(TSig,FSamp,40,0,0,0,1);
%       % Compare with Matlab "Spectrum" Process
%       h = spectrum.welch;
%       CPSD = psd(h,TSig,'FS',FSamp);
%       figure;
%       plot(CPSD)
%       % Real-Valued Sig PSD ... Freq Fold-Over of f = 7300
%       [RSigPSD,RFScale,RSigPow]=SigPSD_V4(real(TSig),FSamp,40,0,0,0,1);
%       % Compare with Matlab "Spectrum" Process
%       h = spectrum.welch;
%       RPSD = psd(h,real(TSig),'FS',FSamp);
%       figure;
%       plot(RPSD)
%
% =========================================================================

function [SigPsd,FScale,SigPow] = SigPSD_V6(SigIn,FSamp,PSDFlr,...
                                              SAve,NAve,SNorm,SPlot)                                          
SigPsd=[];
FScale=[];
SigPow=[];
[NSigs,NSamps]=size(SigIn);
SigType=isreal(SigIn); % SigType=1 -> Real , SigType=0 -> Complex

% Single Column Vector Input ? ... Yes -> Transpose & Proceed
NTFlg=0; % Initialize Transpose Flag
if NSamps==1&&NSigs>=1
    NTFlg=1;
    SigIn=SigIn.';
    [NSigs,NSamps]=size(SigIn);
end

if SigType==0 % Complex
    SigIn=conj(SigIn); % [V6] Addition
end

FDim=2; % Perform Row-Wise FFT [Matrix Rows or Row Vector]
SigFFT=fft(SigIn,[],FDim)/NSamps;
TmpPsd=abs(SigFFT).^2; % Non-FFTshifted PSD:  0 < PSD < FSamp

SigPow=sum(TmpPsd,2);

if(SAve==1) % Apply Row-Wise Averaging/Smoothing to PSD
    NTest = mod(NAve,2); % Check/Adjust NAve to be an Odd Number
    if NTest==0 % Even NAve Input
        NAve=NAve+1; % Adjust NAve to be Odd
    end
    %
    for k=1:NSigs
        CurrPsd=TmpPsd(k,:);
        TmpPsd(k,:)=smooth(CurrPsd,NAve);
    end
end

% Set PSD Floor Level to 'PSDFlr' Below Max[PSD]
if PSDFlr~=0 % Set Floor Value(s) [V5]
    PSDFlr=abs(PSDFlr); % Accept +/-PSDFlr Input
    PSDmax=max(TmpPsd.');
    PSDFlrRat=10^(PSDFlr/10);
    PSDfloor=PSDmax/PSDFlrRat;
    for m=1:NSigs % Threshold / Set PSD Floor Values
        [~,FlrDex]=find(TmpPsd(m,:)<PSDfloor(m));
        TmpPsd(m,FlrDex)=PSDfloor(m);
    end
end

% Normalize PSD(s) if Selected
if SNorm==1
    for m=1:NSigs
        PSDmax=max(TmpPsd(m,:).'); % Update PSDmax
        TmpPsd(m,:)=TmpPsd(m,:)./PSDmax;
    end
end

% Update Output PSD
if NTFlg==1; % Return Output to Column Vector
    SigPsd=TmpPsd.';
else
    SigPsd=TmpPsd;
end

DelFrq=FSamp/NSamps;
FScale=(0:DelFrq:(NSamps-1)*DelFrq); % Output PSD Freq Scale

if SPlot==1 % Produce Output Plots
    disp(' ')
    disp('SigPSD_V6: Generating Plots')
    disp(' ')
    %
    if NSigs==1
        SigPlot=SigIn;
        FFTPlot=SigFFT;
        PsdPlot=SigPsd;
    else % Use 1st SigIn & SigPsd Matrix for Plot Default
        SigPlot=SigIn(1,:);
        FFTPlot=SigFFT(1,:);
        PsdPlot=SigPsd(1,:);
    end
    %
    % Get Current Figure Number
    %
    CurFigHand=get(0,'Children');
    if isempty(CurFigHand)==1
        FigDex=1;
    else
        FigMaxDex=0;
        [NCurFig,~]=size(CurFigHand);
        for m=1:NCurFig
            if CurFigHand(m,:).Number>FigMaxDex
                FigMaxDex=CurFigHand(m,:).Number;
            end
            FigDex=FigMaxDex+1;
        end
    end
    %
    figure(FigDex)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,1) % Top Plot: Time Domain Response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TimScale=linspace(0,NSamps/FSamp,NSamps);
    if SigType==1 % SigIn is Real
        plot(TimScale, SigPlot,'b');
        grid;
        set(gca,'XLim',[0 max(TimScale)]);
        set(gca,'YLim',[-1.1*max(abs(SigPlot)) 1.1*max(abs(SigPlot))]);
        if NSigs==1
            title([{'SigPSD-V6: Real-Valued Input Signal'},...
                {'Time Domain Response'}])
        else % NSigs > 1
            title([{['SigPSD-V6: ',num2str(NSigs),...
                ' Real-Valued Input Signals']},...
                {['Signal #1 Time Domain Response']}])
        end
        xlabel('Time (Sec)')
        ylabel('Amplitude')
    else % SigIn is Complex ... Plot Magnitude as Abs(SigPlot)
        plot(TimScale, abs(SigPlot),'b');
        grid;
        set(gca,'XLim',[0 max(TimScale)]);
        if NSigs==1
            title([{'SigPSD-V6: Complex-Valued Input Signal'},...
                {'Time Domain Magnitude Response'}])
        else % NSigs > 1
            title([{['SigPSD-V6: ',num2str(NSigs),...
                ' Complex-Valued Input Signals']},...
                {['Time Domain Magnitude Response for Signal #1']}])
        end
        xlabel('Time (Sec)')
        ylabel('Abs[s(t)]')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,2) % Middle Plot:  PSD Response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PsdPlotDb=10*log10(PsdPlot);
    
    if SigType==1 % SigIn Real: Plot FFTshift{PSD} for -Fs/2 < f < Fs/2
        FrqPlot=FScale-FSamp/2;
        plot(FrqPlot,fftshift(PsdPlotDb),'b'); % PSD Shift
        xlabel('-FSamp/2 < f < FSamp/2 (Hz)')
    else % SigIn Complex: Plot PSD for 0 < f < Fs
        FrqPlot=[0:DelFrq:(NSamps-1)*DelFrq];
        plot(FrqPlot,PsdPlotDb,'b'); % No PSD Shift
        xlabel('0 < f < FSamp (Hz)')
        ylabel('| FFT |^ 2')
    end
    ylabel('| FFT |^ 2')
    grid
    % Adjust X & Y Axis Limits
    PsdMax=max(PsdPlotDb);
    PsdMin=min(PsdPlotDb);
    YPsdDel=0.10*(PsdMax-PsdMin); % Arb 10%
    set(gca,'YLim',[PsdMin-YPsdDel PsdMax+YPsdDel])
    set(gca,'XLim',[FrqPlot(1) FrqPlot(end)])
    
    % Add Normalized / Un-Normalized Titles
    if SNorm==1 % Normalized
        title('Normalized PSD Response')
    else % Unormalized
        title('Un-Normalized PSD Response')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,3) % Bottom Plot: Angle of FFT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FFTang=angle(FFTPlot);
    if SigType==1 % Real-Valued Input
        plot(FrqPlot, fftshift(unwrap(FFTang)),'b')
        xlabel('-FSamp/2 < f < FSamp/2 (Hz)')
    else
        plot(FrqPlot, unwrap(FFTang),'b')
        xlabel('0 < f < FSamp (Hz)')
    end
    grid
    AngMax=max(unwrap(FFTang));
    AngMin=min(unwrap(FFTang));
    YAngDel=0.10*(AngMax-AngMin); % Arb 10%
    set(gca,'YLim',[AngMin-YAngDel AngMax+YAngDel])
    set(gca,'XLim',[FrqPlot(1) FrqPlot(end)])
    title('Unwrapped FFT Phase Response')
    ylabel('Angle (Rads)')
    
end % SigPSD Function

% ========================== End SigPSD Function ==========================