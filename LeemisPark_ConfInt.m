% =========================================================================
%
%      function [CI]= LeemisPark_ConfInt(MeanIn,NumMC,NPlot)
%
% =========================================================================
%
% Created: 1 Nov 2009 (Dr. Temple)
% Modified: Feb 2010 (Dr. Temple)
%
% Confidence Intervals (CI) Calculated per Leemis & Park
%   Ref: Leemis, L.M. and S.K. Park. Discrete-Event Simulation: A First
%        Course. Prentice Hall, New Jersey, 2006.
%
%   CI: MeanIn +/- 1.96*sqrt[MeanIn(1-MeanIn)/NumMC]
%
%
% Inputs:   MeanIn: Mean Values (Vector)
%           NumMC: # Monte Carlo Trials for Each Mean (Vector)
%           Nplot: Plot Cntl Variable ... 1 -> Plot  &  0 -> No Plot
%
% Input Test Data
%
% EbNoDb=[-2:.5:5]; % Db
% EbNo=10.^(EbNoDb/10); % Ratio Form
% MeanIn=qfunc(sqrt(2*EbNo))
% NumMC=ones(1,length(EbNo))*1000
% [CI]= LeePark_ConfInt(MeanIn,NumMC,1);
%
%
% =========================================================================

function [CI]= LeemisPark_ConfInt(MeanIn,NumMC,NPlot)

[NumMeans]=length(MeanIn);

for k=1: NumMeans
    CI(k) = 1.96*sqrt((MeanIn(k)*(1-MeanIn(k)))/NumMC(k));
end

if NPlot==1
    % Plot MeanIn with CIs
    figure
    hold
    for k=1:3
        errorbar(MeanIn,CI,'-*');
    end
    grid
    set(gca,'YScale','log') % Set Y-Scale to Log10 vs. Linear
    title('Plot of MeanIn Values with CIs')
    ylabel('MeanIn With CIs' )
    xlabel('MeanIn Number')
end

% =========================================================================

