% %% LOAD DATA %% %
params = proj_config();

% Load DataStruct
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

% Load Network and Metrics data for baseline, pre cue, and pre go

%{
Baseline = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('baseline_%s_.mat', params.subj)));
%}
preCue = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_1_%s_.mat', params.subj)));
preGo = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_2_%s_.mat', params.subj)));


% Load Loc data for patient. 
load(params.subjChLoc);
i_pt = strcmp(params.subj, {patient_loc(1).session.subjID});
ptChLoc = patient_loc(1).session(i_pt);
if sum(i_pt) > 1
    ptChLoc = ptChLoc(1);
end

ichan = (ptChLoc.type~=0);
coords = ptChLoc.coords(ichan,:); 
cnames = ptChLoc.names(ichan,:); 

Ntrl = dataStruct.Ntrl;
%Xcue = ones(Ntrl,dataStruct.Nch,5);
%Xgo = ones(Ntrl,dataStruct.Nch);
X = ones(Ntrl,dataStruct.Nch,6);
y = dataStruct.RT;
for i = 1:Ntrl
    X(i,:,1:5) = preCue.Metrics(i).aveCtrl(:,:);
    X(i,:,6) = preGo.Metrics(i).aveCtrl(:,end);   
end

meanX = mean(X);
stdX = std(X);
SEM = stdX/sqrt(Ntrl); 
CI95 = tinv([0.025 0.975], Ntrl-1);
XCI95 = bsxfun(@times, SEM, CI95(:));

startInd = 1;
colors = ['r','g','b','m','y','c','k'];
for i=1:8
    if i <= rem(Ntrl,8)
        endInd = startInd + ceil(Ntrl/8) - 1;
    else
        endInd = startInd + floor(Ntrl/8) - 1;
    end
    figureo
    hold on
    for j = startInd:endInd
        thisMean = squeeze(meanX(1,j,:));
        thisCI = squeeze(XCI95(:,j,:));
        c = colors(mod(j,length(colors))+1);
        plot(1:6, thisMean, 'Color', c, 'DisplayName', num2str(j));
        fill(1:6, thisCI, c);
    end
    hold off
    legend();
    startInd = endInd + 1;
end

%{
Xcue = preCue.Metrics(cleanTrialIdxs(i)).aveCtrl(:,end);
Xgo(i,:) = preGo.Metrics(cleanTrialIdxs(i)).aveCtrl(:,end);
%}