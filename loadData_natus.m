function [dat,Nsamp,fs,uch,chnm] = loadData_natus(varargin)
% function [dat,Nsamp,fs,uch,chnm] = loadData_natus(varargin)
%   Load data recorded on Natus acquisition system in the EMU
%
%   dat = samples x channels (concatonates samples across all session times)
%   Nsamp = number of samples in each session time (so can break out individual sessions)
%   fs = sampling rate (samples/s)
%   uch = channel numbers
%   chnm = channel names
%
%   DR 10/2018

% parameters
ddir = '/Users/drew/Data/CCDT/eeg/'; % data directory
subj = 'HUP152'; % subject
sess = '10Nov17'; % session date
stime = []; % session time (leave empty for all session times)
reref = 'noreref'; % 'reref' or 'noreref'
verbose = 1; 
if nargin
    v2struct(varargin{1});
end

% load params.txt
pth = [ddir subj '/eeg.' reref];
fid = fopen(fullfile(pth, 'params.txt'));
C = textscan(fid,'%s %s');
fclose(fid);
if strcmp(C{1}{1},'samplerate')
    fs = str2double(C{2}{1});
    if verbose, disp(['samplingrate = ' num2str(fs)]); end
else
    error('sampling rate not loaded from params.txt');
end
if strcmp(C{1}{2},'dataformat')
    dform = C{2}{2}(2:end-1);
    if verbose, disp(['dataformat = ' dform]); end
else
    error('data format not loaded from params.txt');
end
if strcmp(C{1}{3},'gain')
    gain = str2double(C{2}{3});
    if verbose, disp(['gain = ' num2str(gain)]); end
else
    error('gain not loaded from params.txt');
end

% load jacksheet.txt
jpth = fopen(fullfile([ddir subj], 'docs/jacksheet.txt'));
J = textscan(fid,'%d %s');
fclose(fid);
uch = J{1}; % channel numbers in jacksheet
chnm = J{2}; % channel names in jacksheet

% parse file names
if ~isempty(stime)
    fnm = dir(fullfile(pth,[subj '_' sess '_' stime '*']));
else
    fnm = dir(fullfile(pth,[subj '_' sess '*']));
end
Nfile = length(fnm);
chnum = []; stimes = [];
for ii = 1:Nfile
    cfnm = fnm(ii).name;
    chnum = [chnum; cfnm(end-2:end)]; % channel number in file name
    stimes = [stimes; cfnm(end-7:end-4)]; % session time in file name
end
uchfnm = unique(chnum,'rows');
uchfnm = str2double(cellstr(uchfnm));
assert(all(uch==uchfnm),'channel numbers in jacksheet and filenames don''t match up');
Nch = length(uch); % number of channels
ustimes = unique(stimes,'rows');
ustimes = str2double(cellstr(ustimes)); 
Nst = length(ustimes); % number of session times

% load data
datc = cell(Nst,1);
for ii = 1:Nfile
    cfnm = fnm(ii).name;
    cch = cfnm(end-2:end); % channel number
    ich = find(uch==str2double(cch));
    cst = cfnm(end-7:end-4); % session time
    ist = find(ustimes==str2double(cst));
    fid = fopen(fullfile(pth, cfnm),'r');
    cdat = fread(fid,inf,dform)*gain;
    fclose(fid);
    for jj = 1:Nst % initialize
        if isempty(datc{jj}) && ist==jj
            datc{jj} = zeros(length(cdat),Nch);
        end
    end
    
    % compile across channels
    datc{ist}(:,ich) = cdat;
end
Nsamp = [];
for ii = 1:Nst-1
    Nsamp = [Nsamp; size(datc{ii},1)]; % samples in each session time
end
dat = cat(1,datc{:}); % concatenate data from multiple session times if more than one
