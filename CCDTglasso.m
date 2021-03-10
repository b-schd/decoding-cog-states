function CCDTglasso
% function CCDTglasso
% Correlate graph metrics with RT performance
%   GN 06/2020

% parameters
verbose = false; %if true, prints a statement for each iteration of GLASSO
pdir = '/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/'; % processed data directory
WstDir = 'WstGLASSO/'; %create this directory to save Wst data
%graphDir = 'graphGLASSO/'; %create this directory to save graph metrics
svnm = 'GLASSO_allF_10secpre_sIall';
p.ddir = '/Users/graceng/OneDrive/CCDTnewDat/CCDT/eeg/'; % data directory
%behavfile = '/Volumes/CCDT/procDataFiles/allSubjBehaviorStruct.mat';


sI = 1; %0 for sI0, 1 for sIall
numRho = 100;
rhoRange = [0.01 0.1];
gamma = 0.5; %0.1 % goodness of fit index for EBIC; must be a value between [0, 1]
maxIt = 1e2;
p.subj = []; % subject (leave empty to batch process all subjects in database)
p.sess = []; % session date (leave empty to chose automatically)
p.stime = []; % session time (leave empty for all session times)
p.rln = 1; % remove line noise? (0 or 1)
p.rrf = 1; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
p.outl = 1; % remove outlier (disconnected) channels? (0 or 1)
p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
p.outlThresh = 5; % outlier threshold (standard deviations)
iev = 1; % event (1 = cue, 2 = go, 3 = response)
win = [-10000 0]; % perievent window (ms)
stsubj = 1; % start subject
%regr = 1; % robust regression? (0 or 1)

% create filepaths for saving data
if not(isfolder([pdir WstDir]))
    mkdir([pdir WstDir])
end
%{
if not(isfolder([pdir graphDir]))
    mkdir([pdir graphDir])
end
%}

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

%{
load(behavfile,{'behav','behav_sIall'});
if sI
    cbehav = behav_sIall;
else
    cbehav = behav;
end
%}

WstStruct = struct;
%{
LTmod = zeros(Nsubj,4);
LTstr = zeros(Nsubj,4);
LTqexp = cell(Nsubj,1);
sP = cell(Nsubj);
sPval = cell(Nsubj);
NFstruct = struct;
SVMstruct = struct;
behavStruct = struct;
%}

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; q.subj = p.subj;
    p.sess = db{isubj,2}; q.sess = p.sess;
    q.ddir = [p.ddir(1:end-4) 'events/'];
    % load data
    [dat,Nsamp,fs,gch, gchlbl] = loadCCDTdata(p);
    Nch = size(dat,2);
    
    % load events
    ccdt = parseCCDTevents(q);
    if ~isempty(Nsamp) && iscell(ccdt)
        nccdt = ccdt{1};
        if sI
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
    
    
    % Compute an effective connectivity matrix for each trial
    allEffConn = zeros(Nch, Nch, Ntrl);
    convergences = zeros(Ntrl,1);
    rhoSet = linspace(rhoRange(1), rhoRange(2), numRho);
    for ii = 1:Ntrl
        allTheta = zeros(numRho, Nch, Nch);
        allEBIC = zeros(numRho,1);
        allConverge = zeros(numRho,1);
        for jj = 1:numRho
            % Input to GraphicalLasso is a Nsamp x Nch matrix
            % Need to pick a rho. Can try 100 values between 0.01 and 0.1,
            % and then use EBIC to pick the optimal value
            % Desired output: theta (dimensions are Nch x Nch)
            % May need to then standardize theta (regularized inverse
            % covariance matrix) to get an effective conn matrix: -1 * omega_i,j
            % / [sqrt(omega_i,i) * sqrt(omega_j,j)]
            S = cov(transpose(squeeze(datwin(:,ii,:))));
            [Theta, ~, logL, converge] = graphicalLassoMatlab(S, rhoSet(jj), verbose, maxIt);
            %need to think about what to do if it doesn't converge -- keep
            %track of ones that don't converge? Try different lambda
            %parameters?
            allTheta(jj,:,:) = Theta;
            allConverge(jj) = converge;
            %a_rho is the number of non-zero off-diagonal elements in Theta(rho)
            a_rho = nnz(Theta - diag(diag(Theta))); %may not work if there is inf or NaN in Theta
            allEBIC(jj) = -2*logL + a_rho*(log(Nsamp) + 4*gamma*log(Nch));
        end
        [~, idxMin] = min(allEBIC);
        bestTheta = squeeze(allTheta(idxMin,:,:));
        allEffConn(:,:,ii) = -1. * bestTheta ./ sqrt(diag(bestTheta)) ./ sqrt(diag(bestTheta)).';
        convergences(ii) = allConverge(idxMin);
    end
    
    Wst = allEffConn;
    WstStruct(isubj).Wst = Wst;
    WstStruct(isubj).convergences = convergences;
    save([pdir WstDir svnm, '.mat'],'WstStruct');
    
    %{
    nodeStr = zeros(Nch, Ntrl);
    for jj = 1:Ntrl
        nodeStr(1:Nch,jj) = sum(Wst(:,:,jj));
    end
    
    NFstr = mean(nodeStr,1);
    zNFstr = zscore(NFstr);
    znodeStr = zscore(nodeStr);

    iF = cbehav(isubj).ifast;
    iS = cbehav(isubj).islow;
    iL = cbehav(isubj).ilate;
    iE = cbehav(isubj).iearly;

    % Graph metrics for each subject per trial
    for i=1:size(Wst,3)
        NFqexp(:,:,i) = getCommunicability(Wst(:,:,i),1,1);
        [Ci(:,i),NFmod(:,i)] = community_louvain(Wst(:,:,i));
    end
    NFqexp_nodal = squeeze(mean(NFqexp,1));

    for i=1:size(NFqexp,3)
       [Ciq(:,i),NFmodq(:,i)] = community_louvain(NFqexp(:,:,i)); 
    end

    % flag trials with bad RTs
    inr = (RT<0 | RT>999); % incorrect responses, logical index
    disp([num2str(sum(inr)) ' trials omitted from analysis']);

    if regr
        for ix=1:length(gch)
            [bqexp(ix,:),statsqexp(ix,:)] = robustfit(NFqexp_nodal(ix,~inr),RT(~inr));
        end
        [bmod,statsmod] = robustfit(NFmod(~inr),RT(~inr));
        [bstr,statsstr] = robustfit(NFstr(~inr),RT(~inr));
        LTmod(isubj,:) = [bmod(2) bmod(2)-1.96*statsmod.se(2) bmod(2)+1.96*statsmod.se(2) statsmod.p(2)];
        LTstr(isubj,:) = [bstr(2) bstr(2)-1.96*statsstr.se(2) bstr(2)+1.96*statsstr.se(2) statsstr.p(2)];
        for ix=1:length(gch)
            LTqexp{isubj}(ix,:) = [bqexp(ix,2) bqexp(ix,2)-1.96*statsqexp(ix).se(2) bqexp(ix,2)+1.96*statsqexp(ix).se(2) statsqexp(ix).p(2)];
        end
        [sPval{isubj},sP{isubj}] = sort(LTqexp{isubj}(:,4), 'ascend');
    end

    % Differential Network Modularity
    try
        [Cidiff,NFmoddiff] = community_louvain(mean(Wst(:,:,iF),3)-mean(Wst(:,:,iS),3));
    catch
        Cidiff = []; NFmoddiff = [];
    end

    try
        [Cidiffq,NFmoddiffq] = community_louvain(mean(NFqexp(:,:,iF),3)-mean(NFqexp(:,:,iS),3));
    catch
        Cidiffq = []; NFmoddiffq = [];
    end

    %SVM
    cN = cell(size(Wst,3),1);
    cN(iF) = {'Fast'};
    cN(iS) = {'Slow'};
    cN(~iF&~iS) = [];

    try
        mdlSVMcv = fitcsvm(NFqexp_nodal((LTqexp{isubj}(:,4)<=0.05),(iF|iS))',cN,'KFold',5);
        SVMnodesUsed = gch(LTqexp{isubj}(:,4)<=0.05);
    catch
        mdlSVMcv = fitcsvm(NFqexp_nodal(sP{isubj}(1),(iF|iS))',cN,'KFold',5);
        SVMnodesUsed = gch(sP{isubj}(1));
    end
    [score_lbl,score_svm] = kfoldPredict(mdlSVMcv);
    [X,Y,T,AUC] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');

    % null comparison
    nullcN = cN(randperm(length(cN)));
    try
        mdlSVMcvnull = fitcsvm(NFqexp_nodal((LTqexp{isubj}(:,4)<=0.05),(iF|iS))',nullcN,'KFold',5);
    catch
        mdlSVMcvnull = fitcsvm(NFqexp_nodal(sP{isubj}(1),(iF|iS))',nullcN,'KFold',5);
    end
    [score_lblnull,score_svmnull] = kfoldPredict(mdlSVMcvnull);
    [Xnull,Ynull,Tnull,AUCnull] = perfcurve(nullcN,score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
    
    NFstruct(isubj).NFqexp_nodal = NFqexp_nodal;
    NFstruct(isubj).NFstr = NFstr;
    NFstruct(isubj).NFmod = NFmod;
    NFstruct(isubj).Ci = Ci;
    NFstruct(isubj).NFmoddiff = NFmoddiff;
    NFstruct(isubj).Cidiff = Cidiff;
    NFstruct(isubj).NFmodq = NFmodq;
    NFstruct(isubj).Ciq = Ciq;
    NFstruct(isubj).NFmoddiffq = NFmoddiffq;
    NFstruct(isubj).Cidiffq = Cidiffq;
    NFstruct(isubj).gch = gch;
    NFstruct(isubj).gchlbl = gchlbl;

    SVMstruct(isubj).mdl = mdlSVMcv;
    SVMstruct(isubj).nodes = SVMnodesUsed;
    SVMstruct(isubj).score_lbl = score_lbl;
    SVMstruct(isubj).score_svm = score_svm;
    SVMstruct(isubj).classNames = cN;
    SVMstruct(isubj).RT = RT(iF|iS);
    SVMstruct(isubj).AUC = AUC;
    SVMstruct(isubj).X = X;
    SVMstruct(isubj).Y = Y;
    SVMstruct(isubj).T = T;
    SVMstruct(isubj).null.mdl = mdlSVMcvnull;
    SVMstruct(isubj).null.score_lbl = score_lblnull;
    SVMstruct(isubj).null.score_svm = score_svmnull;
    SVMstruct(isubj).null.classNames = nullcN;
    SVMstruct(isubj).null.AUC = AUCnull;
    SVMstruct(isubj).null.X = Xnull;
    SVMstruct(isubj).null.Y = Ynull;
    SVMstruct(isubj).null.T = Tnull;

    behavStruct(isubj).subj = db{isubj,1};
    behavStruct(isubj).sess = db{isubj,2};
    behavStruct(isubj).vRT = RT;
    behavStruct(isubj).vDT = DT;
    behavStruct(isubj).ifast = iF;
    behavStruct(isubj).islow = iS;
    behavStruct(isubj).iearly = iE;
    behavStruct(isubj).ilate = iL;
    behavStruct(isubj).error = inr;

    clear NFq* NFstr NFmod* Ci* bmod bqexp bstr stats* mdlSVMcv cN score* X Y T AUC Xnull Ynull Tnull AUCnull SVMnodesUsed
    save([pdir graphDir svnm,'.mat'],'LT*','sP*','NFstruct', 'SVMstruct', 'behavStruct')
    %}
end

end