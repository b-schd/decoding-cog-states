% parameters
cleanTrials = true; %if true, remove all trials with negative reaction times

% %% LOAD DATA %% %

params = proj_config();

% Load DataStruct
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

% Load Network and Metrics data for baseline, pre cue, and pre go
Baseline = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('baseline_%s_%s.mat', params.subj, dataStruct.sessName)));

preCue = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_1_%s_%s.mat', params.subj, dataStruct.sessName)));

preGo = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_2_%s_%s.mat', params.subj, dataStruct.sessName)));


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

% Regression using average controllability in only the last 500ms before the cue
if cleanTrials
    Ntrl = sum(dataStruct.RT >= 0);
    cleanTrialIdxs = find(dataStruct.RT >= 0);
    X = ones(Ntrl,dataStruct.Nch);
    y = dataStruct.RT(cleanTrialIdxs);
    for i = 1:Ntrl
        X(i,:) = preCue.Metrics(cleanTrialIdxs(i)).aveCtrl(:,5);
    end
else
    Ntrl = dataStruct.Ntrl;
    X = ones(Ntrl,dataStruct.Nch);
    y = dataStruct.RT;
    for i = 1:Ntrl
        X(i,:) = preCue.Metrics(i).aveCtrl(:,5);
    end
end

%
% basic multiple linear regression
% mdl = fitlm(X,y);

% FEATURE SELECTION
% stepwise
    %alternative: sequentialfs does sequential selection of features and
    %cross validation

    %TODO: use 5-fold cross-validation to do model selection (pick a 
    %subset of the channels). Find the one that has the best performance,
    %then train that collection of channels on the whole dataset and record
    %adj_Rsquared

% step_mdl = stepwiselm(X,y, 'PEnter',0.06);


% Lasso regularization
%{
lasso_k_mdl = fitrlinear(X, y, 'KFold', 5, 'Lambda',10.^(-(10:-2:2)), ...
    'Learner', 'leastsquares', 'Regularization', 'lasso');
%}

[B,FitInfo] = lasso(X, y, 'CV', 5); %'Lambda',10.^(-(10:-2:2))
    % by default, lasso checks a geometric sequence of Lambdas
idxLambdaMinMSE = FitInfo.IndexMinMSE;
idxLambda1SE = FitInfo.Index1SE;
lassoPlot(B, FitInfo, 'PlotType', 'CV');
legend('show');
%why do we want to use the largest lambda within 1 SD of MSE, instead of using lambda with the lowest MSE?
coef = B(:,idxLambda1SE); 
coef0 = FitInfo.Intercept(idxLambda1SE);
yhat = X*coef + coef0; %fit the best model on the whole dataset
figure(2);
hold on
scatter(y, yhat)
plot(y, yhat)
xlabel('Actual Reaction Times')
ylabel('Predicted Reaction Times')
hold off
yresid = y - yhat;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(B)-1);
% SAD FINDING: if the best cross-validated prediction error is achieved
% when all coefficients are 0, then you should suspect that no linear
% combination of any subset of the regressors may be useful for predicting
% the outcomes


% SVM -- check Vivek's scripts for reference

% PLS regression (uses PLSR)

%PRIORITY:
%bucket fast and slow trials? or fast, moderate, and slow trials and do
%logistic regression

%look at all 2500ms precue
%look at preGo

%repeat analysis with strength and with modal controllability