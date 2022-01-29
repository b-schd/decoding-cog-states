function [ccdt,fs] = parseCCDTevents(varargin)
% function [ccdt,fs] = parseCCDTevents(varargin)
%   Parse data from CCDT events file. Returns matrix:
%   ccdt = trial x [cue indices, go indices, response indices, DT (ms), RT (ms)]
%   fs = sampling rate
%
%   DR 10/2018

% parameters
ddir = '/Volumes/CCDT/events/'; % data directory
subj = 'HUP069'; % subject
sess = '04Oct17'; % session date
stime = []; % session time
if nargin
    v2struct(varargin{1});
end

% load events
startingPath = pwd;
cd(ddir);
load([subj '_events.mat'],'-mat');
eval(['evdata = events;']); clear events % built-in function called events
N = length(evdata); % number of events

% find sampling rate from params.txt (only needed if re-computing DT/RT)
cd ..
cd(['./eeg/' subj '/eeg.noreref']);
fid = fopen('params.txt');
C = textscan(fid,'%s %s');
fclose(fid);
if strcmp(C{1}{1},'samplerate')
    fs = str2double(C{2}{1});
    disp(['samplingrate = ' num2str(fs)]);
else
    error('sampling rate not loaded from params.txt');
end

% find number of sessions and trials
titr = [];
for ii = 1:N
    cfnm = evdata(ii).eegfile;
    if ~strcmp(sess,cfnm(end-11:end-5))
        continue;
    end
    titr = [titr; str2double(cfnm(end-3:end)) evdata(ii).trial];
end
ustimes = unique(titr(:,1));
Nst = length(ustimes); % number of session times
ccdt = cell(Nst,1);
for ii = 1:Nst
    ind = find(titr(:,1)==ustimes(ii));
    ccdt{ii} = zeros(max(titr(ind,2)),5); % initialize output (one cell per session time)
end

% check that events are associated with specified session
for ii = 1:N
    cfnm = evdata(ii).eegfile;
    cst = cfnm(end-3:end); % current session time
    ist = find(ustimes==str2double(cst));
    ind = strfind(cfnm,'/');
    if ~isempty(stime)
        cfnm = cfnm(ind(end)+1:end);
        dfnm = [subj '_' sess '_' stime];
    else
        cfnm = cfnm(ind(end)+1:end-5);
        dfnm = [subj '_' sess];
    end
    if ~strcmp(cfnm,dfnm)
%         warning(['skipped event due to subject/session identifier: ' cfnm]);
        continue;
    end
    
    % extract timestamps (cue, go, response), DT, and RT for this session
    ctrl = evdata(ii).trial;
    if strcmp(evdata(ii).type,'FIX_START') % cue indices
        ccdt{ist}(ctrl,1) = evdata(ii).eegoffset;
    elseif strcmp(evdata(ii).type,'CC') % go indices
        ccdt{ist}(ctrl,2) = evdata(ii).eegoffset;
    elseif strcmp(evdata(ii).type,'RESPONSE') % response indices
        ccdt{ist}(ctrl,3) = evdata(ii).eegoffset;     
    end
    ccdt{ist}(ctrl,4) = (ccdt{ist}(ctrl,2)-ccdt{ist}(ctrl,1))/fs*1000; % delay time (ms; same as evdata(ii).delay)
    ccdt{ist}(ctrl,5) = (ccdt{ist}(ctrl,3)-ccdt{ist}(ctrl,2))/fs*1000; % reaction time (ms; same as evdata(ii).rt)
end
if ~isempty(stime)
    ist = find(ustimes==str2double(stime));
    ccdt = ccdt{ist};
elseif Nst==1
    ccdt = ccdt{1};
end
if strcmp(subj,'HUP069')
    ccdt(42:45,:) = []; % recording blanked out for these 4 trials
end
cd(startingPath);
