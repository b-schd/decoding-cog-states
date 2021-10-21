function [dat,Nsamp,fs,chn,chnm] = loadCCDTdata(varargin)
% function [dat,Nsamp,fs,chn,chnm] = loadCCDTdata(varargin)
%   Load data recorded on Natus acquisition system in the EMU during CCDT
%   task. 
%
%   DR 03/2020

% parameters
ddir = '/Volumes/CCDT/CCDT/eeg/'; % data directory
subj = 'HUP142'; % subject
sess = '30Jun17'; % session date
stime = []; % session time (leave empty for all session times)
rln = 1; % remove line noise? (0 or 1)
rrf = 2; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
outl = 1; % remove outlier channels? (0 or 1)
bpfilt = [.5, 150]; % band pass filter cuttoffs (Hz) (leave empty if no filter)
outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
outlThresh = 5; % outlier threshold (standard deviations)
if nargin
    v2struct(varargin{1});
end

% load data
p.ddir = ddir; p.subj = subj; p.sess = sess; p.stime = stime; p.verbose = 0;
[dat,Nsamp,fs,chn,chnm] = loadData_natus(p);

% remove non-sEEG (DC, EKG, EOG, EEG) channels
inseeg = find(strncmp(chnm,'DC',2)|strncmp(chnm,'EKG',3)|strcmp(chnm,'LOC')|strcmp(chnm,'ROC')|...
              strncmp(chnm,'C',1)|strncmp(chnm,'F',1)|strncmp(chnm,'P',1)|...
              strncmp(chnm,'O',1)|strncmp(chnm,'T',1)|strncmp(chnm,'FP',2));
if ~isempty(inseeg)
    dat(:,inseeg) = [];
    chn(inseeg) = [];
    disp([num2str(length(inseeg)) ' non-sEEG channels removed']);
    disp(chnm(inseeg)');
    chnm(inseeg) = [];
end

% remove null channels (no signal)
dat = dat - ones(size(dat,1),1)*median(dat,1); % remove dc offset
sigmaN = median(abs(dat)/.6745); % standard deviation of background noise
inull = find(sigmaN==0);
if ~isempty(inull)
    dat(:,inull) = [];
    chn(inull) = [];
    chnm(inull) = [];
    disp([num2str(length(inull)) ' null channels removed']);
end

Nch = size(dat,2);

% remove line noise
if rln
    for ii = 1:Nch
        ln = mtmlinenoise(dat(:,ii),3,fs,fs,[60 120 180]);
        dat(:,ii) = dat(:,ii) - ln;
    end
end

% remove outlier channels
if outl
    if strcmp(outlMetric,'rms')
        s = sqrt(median(dat.^2)); % root-median-square amplitude
    elseif strcmp(outlMetric,'powGamma')
        s = bandpower(dat,fs,[80 200]); % high-gamma power (HFA activity)
    elseif strcmp(outlMetric,'kurtosis')
        s = kurtosis(dat); % kurtosis
    end
    zs = 0.6745*(s-median(s))/mad(s,1); % modified zscore of outlier statistic (Iglewicz and Hoaglin 1993)
    ioutlier = find(abs(zs) > outlThresh); % Iglewicz and Hoaglin recommended 3.5 standard deviations
    if ~isempty(ioutlier)
        dat(:,ioutlier) = [];
        chn(ioutlier) = [];
        chnm(ioutlier) = [];
        disp([num2str(length(ioutlier)) ' outlier channels removed']);
    end
end

% common average re-reference
if rrf==1
    dat = dat - mean(dat,2)*ones(1,size(dat,2));
end

% bandpass filter
if bpfilt
    [b,a] = butter(2, bpfilt/fs); 
    dat = filtfilt(b, a, dat);
end

% monopolar, bipolar, laplacian montage
if rrf>=2 % multipolar
    charnm = char(chnm);
    if size(charnm,2)==5
        [leads,~,ilead] = unique(charnm(:,1:3),'rows');
    else
        [leads,~,ilead] = unique(charnm(:,1:2),'rows');
    end
    datrr = zeros(size(dat));
    ind = 1; nchn = []; nchnm = [];
    for ii = 1:size(leads,1)
        ich = find(ilead==ii);
        if rrf==2 && length(ich)>1 % bipolar re-reference: i - i+1
            Ncbp = length(ich)-1;
            datrr(:,ind:ind+Ncbp-1) = dat(:,ich(1:end-1)) - dat(:,ich(2:end));
            ind = ind + Ncbp;
            nchn = [nchn; (chn(ich(1:end-1))+0.5)];
            nchnm = [nchnm; chnm(ich(1:end-1)) chnm(ich(2:end))]; % technically halfway between
        elseif rrf==3 && length(ich)>2 % laplacian re-reference: i - (i-1 + i+1)/2 (see Li et al 2018 NeuroImage)
            Ncbp = length(ich);
            datrr(:,ind) = dat(:,ich(1)) - dat(:,ich(2));
            datrr(:,ind+1:ind+Ncbp-2) = dat(:,ich(2:end-1)) - (dat(:,ich(1:end-2)) + dat(:,ich(3:end)))/2;
            datrr(:,ind+Ncbp-1) = dat(:,ich(end-1)) - dat(:,ich(end));
            ind = ind + Ncbp;
            nchn = [nchn; chn(ich)];
            nchnm = [nchnm; chnm(ich)];
        end
    end
    datrr(:,ind:end) = [];
    dat = datrr;
    chn = nchn;
    chnm = nchnm;
end
