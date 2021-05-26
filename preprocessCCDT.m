function [dataStruct, Nsubj] = preprocessCCDT(ddir, subjChLoc, subj, pdir)
%{
- ddir: data directory ('/Users/graceng/OneDrive/CCDTnewDat/CCDT/eeg/')
- subjChLoc: directory for file with information about channel location
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/'patient_loc.mat')
- subj: enter string name of the subject (eg, 'HUP069') or enter [] to
    process all subjects
- pdir (optional): processed data directory
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/processData/')
%}


% parameters
p.ddir = ddir;
p.subjChLoc = subjChLoc; 
if nargin > 3
    p.pdir = pdir;
else
    p.pdir = [];
end

p.sI = 1; %0 for sI0, 1 for sIall
p.subj = subj;
p.sess = []; % session date (leave empty to chose automatically)
p.stime = []; % session time (leave empty for all session times)
p.rln = 1; % remove line noise? (0 or 1)
p.rrf = 1; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
p.outl = 1; % remove outlier (disconnected) channels? (0 or 1)
p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
p.outlThresh = 5; % outlier threshold (standard deviations)
iev = 1; % event (1 = cue, 2 = go, 3 = response)
win = [-10000 0]; % perievent window (ms) %[-500 0]
stsubj = 1; % start subject

% check filepaths for saving data
if ~isempty(p.pdir)
    if ~isfolder(p.pdir)
        error('pdir does not exist.');
    end
end

% load database
db = CCDTdatabase;
if ~isempty(p.subj)
    ind = cellfun(@(x) strcmp(x,p.subj),db(:,1)); % find specified subj
    if ~isempty(p.sess)
        ind2 = cellfun(@(x) strcmp(x,p.sess),db(:,2)); % find specified session
        ind = ind & ind2;
    end
    ind = find(ind);
    if isempty(ind), error('subj/sess not found in database'); end
    db = db(ind(1),:);
end

% load data
Nsubj = size(db,1);
load(p.subjChLoc);

dataStruct = struct;
for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; q.subj = p.subj;
    p.sess = db{isubj,2}; q.sess = p.sess;
    q.ddir = [p.ddir(1:end-4) 'events/'];
    
    % load data for this subject
    [dat,Nsamp,fs,gch, gchlbl] = loadCCDTdata(p);
    gchlbl(patient_loc(1).session(isubj).type==0)=[];
    gch(patient_loc(1).session(isubj).type==0)=[];
    dat(:,patient_loc(1).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    % load events for this subject
    ccdt = parseCCDTevents(q);
    if ~isempty(Nsamp) && iscell(ccdt)
        nccdt = ccdt{1};
        if p.sI
            disp('Concatenating within-day multisession')
            for ii = 2:length(ccdt)
                cccdt = ccdt{ii};
                cccdt(:,1:3) = cccdt(:,1:3) + sum(Nsamp(1:ii-1))*ones(size(cccdt,1),3); % event index accounting for concatenated data across session times
                nccdt = [nccdt; cccdt];
            end
        end
        ccdt = nccdt;
    end
    Ntrl = size(ccdt,1);
    CT = ccdt(:,1)/fs; % cue time (s)
    DT = ccdt(:,4); % delay time (s)
    RT = ccdt(:,5); % reaction time (s)

    
    % window indices
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin);
    indmat = ccdt(:,iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin;

    % window data
    datwin = zeros(Nch,Ntrl,Nsamp);
    for ii = 1:Nch
        cdat = dat(:,ii);
        datwin(ii,:,:) = cdat(indmat);
    end
    
    dataStruct(isubj).datwin = datwin;
    dataStruct(isubj).Ntrl = Ntrl;
    dataStruct(isubj).Nch = Nch;
    dataStruct(isubj).Nsamp = Nsamp;
    dataStruct(isubj).CT = CT;
    dataStruct(isubj).DT = DT;
    dataStruct(isubj).RT = RT;
    dataStruct(isubj).gchlbl = gchlbl;
    dataStruct(isubj).gch = gch;
    dataStruct(isubj).subjName = db(isubj);
end

if ~isempty(p.pdir)
    save([p.pdir, 'processedData.mat'], 'dataStruct', 'Nsubj');
end

end