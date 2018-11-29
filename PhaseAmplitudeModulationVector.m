% following the definition in Canolty et al. 2006
% code by Paola Malerba - finalized on Nov 29th 2018

% this function considers two different (or same) signals - one will provide the phase and one the amplitude - and finds the phase-amplitude modulation vector
% within the signal, this function cuts only times that are of interest (in my case, during slow oscillations) 
% and discards the parts of signal that are not of interest before computing the modulation vector but after finding phase and amplitude 
% (so no edge effects are inserted)

% INPUT
% phaseSig --> signal that will give the phase (i.e. pre-filtered in the range of interest)
% this is a 1-d array of voltage in time (time binning on phaseSig and ampSig has to be the same)

% ampSig --> signal that will give the amplitude (i.e. pre-filtered in the range of interest)
% this is a 1-d array of voltage in time 

% freq --> value in Hz of the sampling frequency (i.e. inverse of time binning length)

% SO --> structure with 2 fields :
% SO.starts
% SO.ends
% the two fields mark the locations in time where the times of interest begin and end 
% if you want to use the whole signal rather than only selected time intervals, just make sure SO.starts = 1 and SO.ends = length(phaseSig);

% OUTPUT
% phasePDF is the probability distribution function of the phases found (make sure it looks flat)
% it has two fields
% phasePDF.y --> the probability axis
% phasePDF.x --> the phase value axis
% plot(phasePDF.x,phasePDF.y) will show you if you are sampling phase uniformily (flat) or if you have a bias (not good)


% amplitude PDF analogue to phasePDF but for the amplitude values
% again plotting it shows if you have a biased dataset on the amplitude component (you might have to look for a different measure if you do)

% prob dist function of z-signal values
% where the z-signal is the composite complex signal
% it has 3 fields:
% zPDF.xphi values of phases (one axis)
% zPDF.xamp values of amplitudes (one axis)
% zPDF.y  2Dimensional density of signal z over xphi and xamp values 

% ModVect
% ModVect.raw is the length of the average value of complex signal z
% measures phase-amplitude interaction but is not controlled for stochasitcity in the signals 

% to establish if coupling between amp and phase is driven by Amp OR Phase independently or by their actual mutual information
% make surrogate data: from Canolty et al 2006 supp Info
% and get
% ModVect.normalized --> after bootstrapping (= shuffling many many times)

%% NOTE: there are a few hardcoded values that you might want to change for yourself
% phstep = pi/18;  % increase the denominator to refine bin size. too much refinement will lead to very noisy distributions and poor PLI
% surrN = 1000; % number of shuffles you want to compute your normalized ModVect. No need to make this huge start at 200 and grow 


--------------------------------------------------------------------------------------


function [phasePDF,ampPDF,zPDF,ModVect]...
    = PhaseAmplitudeModulationVector...
    (phaseSig,ampSig,freq,SO)

phstep = pi/18;  % increase the denominator to refine bin size. too much refinement will lead to very noisy distributions and poor PLI
phase_x = -pi:phstep:pi;

HphaseSig = hilbert(phaseSig);
HampSig = hilbert(ampSig);

Amp_t = sqrt(HampSig.*conj(HampSig));
Ph_t = angle(HphaseSig); 


% only use portions of the signal which have detected SOs in them
SOtime = zeros(size(Amp_t));
for idso = 1:length(SO.starts)
    % SO.starts is in time-stamps, not time
    % SO.ends is in time-stamps, not time
    xso = SO.starts(idso):SO.ends(idso);
    SOtime(xso) = 1;
end
AmpCut = Amp_t(SOtime>0);
PhCut = Ph_t(SOtime>0);

%outputs
n = hist(PhCut,phase_x);
phasePDF.y = n/sum(n);
phasePDF.x = phase_x;
[n,xout] = hist(AmpCut,100);
ampPDF.y = n/sum(n);
ampPDF.x = xout;

% build interaction signal to study mutual info between amp signal and
% phase signal
Zsig = AmpCut.*exp(1i*PhCut);
zPDF.xphi = phase_x;
zPDF.xamp = xout;
z2D = zeros(length(xout),length(phase_x));
for idA = 1:(length(zPDF.xamp)-1)
    if idA+1<length(zPDF.xamp)
        xA = find(and(AmpCut>=zPDF.xamp(idA), AmpCut<zPDF.xamp(idA+1)));    
    else
        xA = find(and(AmpCut>=zPDF.xamp(idA), AmpCut<=zPDF.xamp(idA+1)));
    end
    zA = Zsig(xA);
    p1 = hist(zA,phase_x);
    z2D(idA,:) = p1;
end
z2D = z2D/sum(z2D(:)); % becomes a density
zPDF.y = z2D;

ModZ = mean(Zsig); % to be compared with surrogate data 
ModVect.raw = ModZ;
% to establish if coupling between amp and phase is driven by Amp OR Phase
% independently or by their actual mutual information

% make surrogate data: from Canolty et al 2006 supp Info
numpts = length(AmpCut);
surrN = 1000; %200;
minSkip = freq;
maxSkip = numpts-freq;
skip = ceil(numpts.*rand(surrN*2,1));
skip(skip>maxSkip) = []; %remove them
skip(skip<minSkip) = [];
skip = skip(1:surrN);
surr_Mod = zeros(surrN,1);
for sID = 1:surrN
    surrAmp = [AmpCut(skip(sID):end) AmpCut(1:skip(sID)-1)];
    zsur = mean(surrAmp.*exp(1i*PhCut));
    surr_Mod(sID) = sqrt(zsur.*conj(zsur));
end
% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox 
[surrMod_MU,surrMod_STD]=normfit(surr_Mod); 
% normalize length using surrogate data (z-score) 
ModNorm_length=(abs(ModZ)-surrMod_MU)/surrMod_STD; 
ModNorm_phase= angle(ModZ); 
ModVect.normalized = ModNorm_length*exp(1i*ModNorm_phase); 

