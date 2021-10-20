function [dataStruct, Nsubj] = preprocessCCDT(ddir, subjChLoc, subj, ...
    params, pdir)
%{
- ddir: data directory ('/Users/graceng/OneDrive/CCDTnewDat/CCDT/eeg/')
- subjChLoc: directory for file with information about channel location
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/patient_loc.mat')
- subj: enter string name of the subject (eg, 'HUP069') or enter [] to
    process all subjects
- params: struct including other parameters such as iev, winLength
- pdir (optional): processed data directory
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/processData/')
%}


% parameters
p.ddir = ddir;
p.subjChLoc = subjChLoc; 
p.subj = subj;
if nargin > 3
    S_out = p;
    fields = fieldnames(params);
    for curr_field = fields(:)'  
      S_out = setfield(S_out,curr_field{:},getfield(params,curr_field{:}));
    end
    p = S_out;
else
    p.stsubj = 1; %index of first subject if processing all subjects
    p.sI = 1;
    p.sess = []; % session date (leave empty to chose automatically)
    p.stime = []; % session time (leave empty for all session times)
    p.reref = 'noreref';
    p.rln = 1; % remove line noise? (0 or 1)
    p.rrf = 1; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
    p.outl = 1; % remove outlier (disconnected) channels? (0 or 1)
    p.bpfilt = [.5, 150]; % band pass filter cuttoffs (Hz) (leave empty if no filter)
    p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
    p.outlThresh = 5; % outlier threshold (standard deviations)
    p.baselineLength = 10000;
    p.iev = 1; % event of interest (1 = cue, 2 = go, 3 = response)
    p.winLength = 500; % time window to analyze prior to event of interest

end
if nargin > 4
    p.pdir = pdir;
else
    p.pdir = [];
end

win = [p.winLength*-1 0]; % perievent window (ms) %[-500 0] or [-10000 0];

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
    p.stsubj = 1;
end

% load data
Nsubj = size(db,1);
load(p.subjChLoc);

dataStruct = struct;
for isubj = p.stsubj:Nsubj
    if isempty(p.subj), ind = isubj; end
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; q.subj = p.subj;
    p.sess = db{isubj,2}; q.sess = p.sess;
    q.ddir = [p.ddir(1:end-4) 'events/'];
    
    % load data for this subject
    [dat,Nsamp,fs,gch, gchlbl] = loadCCDTdata(p);
    % remove bad channels
    gchlbl(patient_loc(1).session(ind).type==0)=[]; %note: patient_loc(1) pulls channel info from the subject's first session
    gch(patient_loc(1).session(ind).type==0)=[];
    dat(:,patient_loc(1).session(ind).type==0) = [];
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

    % get indices for baseline time window
    bwin = round(p.baselineLength*-1/1000*fs):0;
    bNsamp = length(bwin);
    bIndMat = ccdt(1,1)*ones(bNsamp,1) + bwin'; % has dimensions samples x 1

    % get data for baseline time window
    bDatwin = zeros(Nch,bNsamp);
    for ii = 1:Nch
        cdat = dat(:,ii);
        bDatwin(ii,:) = cdat(bIndMat);
    end
    
    % get indices for time window of interest
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin);
    indmat = ccdt(:,p.iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin; 
    % indmat has dimensions Ntrl x Nsamp (samples_per_session). It gives
    % the indices for all samples in the window of time prior to the event
    % denoted by iev

    % get data for time window of interest
    datwin = zeros(Nch,Ntrl,Nsamp);
    for ii = 1:Nch
        cdat = dat(:,ii); % dat has dimensions samples x channels, with all samples concatenated
        datwin(ii,:,:) = cdat(indmat); % looks at all data in the pre-trial period
    end
    
    dataStruct(isubj).bDatwin = bDatwin;
    dataStruct(isubj).datwin = datwin;
    dataStruct(isubj).Ntrl = Ntrl;
    dataStruct(isubj).Nch = Nch;
    dataStruct(isubj).Nsamp = Nsamp;
    dataStruct(isubj).fs = fs; 
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