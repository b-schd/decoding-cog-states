% parameters
cleanTrials = false; %if true, remove all trials with negative reaction times
freqBand = 'High Gamma'; %'none' for no band pass
freqBandMap = containers.Map({'none', 'Theta/Alpha', 'Beta', 'Low Gamma', ...
    'High Gamma'}, {[0 150], [3 12], [14 30], [35 55], [70 150]}); 
            %{} if data is not split up according to frequency bands
plrty = 1; %polarity can be 1 or 2 for channel location file
k = 10; %number of folds for k-fold cross validation
n_rep = 10; %number of times to repeat k-fold cross validation
nboot = 1000;

%%% LOAD DATA %%%
params = proj_config();

% Load DataStruct
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

% Load Network and Metrics data for baseline, pre cue, and pre go
bands = freqBandMap(freqBand);
preCue = load(fullfile(params.ndir, params.subj, 'networks', ...
    sprintf('%s_iev-1_bp-%d-%d.mat', params.subj, bands(1), bands(2))));
preGo = load(fullfile(params.ndir, params.subj, 'networks', ...
    sprintf('%s_iev-2_bp-%d-%d.mat', params.subj, bands(1), bands(2))));

% Load SEEG location data 
load(params.subjChLoc);
db = CCDTdatabase;
subjInd = cellfun(@(x) strcmp(x,params.subj),db(:,1)); % find specified subj
if ~isempty(params.sess)
    subjInd2 = cellfun(@(x) strcmp(x,params.sess),db(:,2)); % find specified session
    subjInd = subjInd & subjInd2;
end
subjInd = find(subjInd);
if isempty(subjInd), error('subj/sess not found in database'); end
locInfo = patient_loc(plrty).session(subjInd);

% Set up average controllability matrices
if cleanTrials
    trialIdxs = find(dataStruct.RT >= 0);
else
    trialIdxs = 1:dataStruct.Ntrl;
end
Ntrl = length(trialIdxs);
Nch = size(preCue.Metrics(1).aveCtrl,1);
precueControl = ones(Ntrl,Nch,size(preCue.Metrics(1).aveCtrl,2));
pregoControl = ones(Ntrl,Nch,size(preGo.Metrics(1).aveCtrl,2));
for i = 1:Ntrl
    precueControl(i,:,:) = preCue.Metrics(trialIdxs(i)).aveCtrl;
    pregoControl(i,:,:) = preGo.Metrics(trialIdxs(i)).aveCtrl;
end
allRT = dataStruct.RT(trialIdxs);

randOrder = randperm(Ntrl);
Xcue = precueControl(randOrder,:,end);
Xgo = pregoControl(randOrder,:,end);
y = allRT(randOrder);
ynull = y(randOrder);
y_binned = (y < median(y)); % split RTs into fast half vs slow half
ynull_binned = (ynull < median(y));

%%% ANALYSES COMPARING EXPERIMENTAL DATA TO NULL MODELS %%%

allX = {Xcue, Xgo};
allY = {y, ynull};
allYbinned = {y_binned, ynull_binned};
allStats = {struct, struct}; %first struct stores stats about preCue data, 
                               % second struct stores stats about preGo data

% Multiple linear regression
allRMSE = {{[], []}, {[], []}};
for m = 1:length(allX)
    thisX = allX{m};
    for n = 1:length(allY)
        thisY = allY{n};
        rmse = zeros(k*n_rep,1);
        for i = 1:n_rep
            CVMdl = fitrlinear(thisX, thisY, 'KFold', k);
            losses = kfoldLoss(CVMdl, 'Mode', 'individual');
            rmse((i-1)*k+1:i*k) = sqrt(losses);
        end
        allRMSE{m}{n} = rmse;
    end
end
for i = 1:length(allStats)
    [~,p,~,~] = ttest2(allRMSE{i}{1}, allRMSE{i}{2});
    allStats{i}.MLR.expRMSE = mean(allRMSE{i}{1});
    allStats{i}.MLR.nullRMSE = mean(allRMSE{i}{2});
    allStats{i}.MLR.p = p;
end


% Logistic regression?

% SVMs
allAUC = {{[], []}, {[], []}};
for m = 1:length(allX)
    thisX = allX{m};
    for n = 1:length(allY)
        thisYbinned = allYbinned{n};
        mdlSVM = fitclinear(thisX, thisYbinned, 'Learner', 'svm', 'KFold', k);
        [label,scores] = kfoldPredict(mdlSVM);
        [~,~,~,auc] = perfcurve(thisYbinned, scores(:,1), 1, 'NBoot', ...
            nboot, 'boottype', 'cper');
        aucSE = (auc(3) - auc(1)) / 1.96;
        allAUC{m}{n} = [auc(1), aucSE]; %stores [mean AUC, AUC standard error]
    end
end
for i = 1:length(allStats)
    x1 = allAUC{i}{1}(1); x2 = allAUC{i}{2}(1);
    s1 = allAUC{i}{1}(2); s2 = allAUC{i}{2}(2);
    tstat = (x1 - x2) / sqrt(s1^2/Ntrl + s2^2/Ntrl); %Welch's formula for 
                                % unpaired t test with unequal variances
    df = (s1^2/Ntrl + s2^2/Ntrl)^2 / ((s1^2/Ntrl)^2/(Ntrl-1) + ...
        (s2^2/Ntrl)^2/(Ntrl-1)); %calculate df using Welch's formula
    p = 2*tcdf(abs(tstat),df,'upper'); %calculate 2-tailed p-value
    allStats{i}.SVM.expAUC = x1;
    allStats{i}.SVM.nullAUC = x2;
    allStats{i}.SVM.p = p;
end


%%%%%%%%%%%%%%%%%%%%% OLD CODE -- needs cleaning %%%%%%%%%%%%%%%%%%%%%%%

%{
% FEATURE SELECTION FOR LINEAR REGRESSION
% stepwise
    %alternative: sequentialfs does sequential selection of features and
    %cross validation

    %TODO: use k-fold cross-validation to do model selection (pick a 
    %subset of the channels). Find the one that has the best performance,
    %then train that collection of channels on the whole dataset and record
    %adj_Rsquared

% step_mdl = stepwiselm(X,y, 'PEnter',0.06);


% Lasso regularization for linear regression

lasso_k_mdl = fitrlinear(X, y, 'KFold', k, 'Lambda',10.^(-(10:-2:2)), ...
    'Learner', 'leastsquares', 'Regularization', 'lasso');

[B,FitInfo] = lasso(X, y, 'CV', k); %'Lambda',10.^(-(10:-2:2))
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



% basic logistic regression
[b1, dev1, stats1] = glmfit(X, y_binned, 'binomial', 'link', 'logit');


% FEATURE SELECTION FOR LOGISTIC REGRESSION
% compare errors using deviance vs. AUC?
% try LASSO and stepwise

% Stepwise feature selection for logistic regression
% edit this for cross validation
startInd = 1;
stepwise_models = {};
stepwise_null_models = {};
for i = 1:k
    if i <= rem(Ntrl,k)
        endInd = startInd + ceil(Ntrl/k) - 1;
    else
        endInd = startInd + floor(Ntrl/k) - 1;
    end
    %Xval = X(startInd:endInd,:); 
    %yval_binned = y_binned(startInd:endInd);
    %yval_null_binned = ynull_binned(startInd:endInd);
    trainInds = setdiff(1:Ntrl,startInd:endInd);
    Xtrain = X(trainInds,:);
    ytrain_binned = y_binned(trainInds);
    ytrain_null_binned = ynull_binned(trainInds);
    mdl = stepwiseglm(Xtrain, ytrain_binned, 'linear', 'Distribution', ...
        'binomial', 'PEnter',0.06);
    null_mdl = stepwiseglm(Xtrain, ytrain_null_binned, 'linear', ...
        'Distribution', 'binomial', 'PEnter',0.06);
    stepwise_models{i} = mdl;
    stepwise_null_models{i} = null_mdl;
    startInd = endInd + 1;
end

% Lasso regularization for logistic regression
[B, FitInfo] = lassoglm(X, y_binned, 'binomial', 'NumLambda', 25, 'CV', k);
lassoPlot(B, FitInfo, 'PlotType', 'CV');
legend('show', 'Location', 'best');
idxLambda1SE = FitInfo.Index1SE;
coef = B(:,idxLambda1SE); 
coef0 = FitInfo.Intercept(idxLambda1SE);
bestB = [coef0;coef];
yhat = glmval(bestB, X, 'logit');
histogram(y_binned - yhat)
title('Residuals from lasso logistic regression')


%{
function pred = classf(X, y_binned)
mdl = fitglm(X, y_binned, modelspec,'Distribution','binomial');
yfit = predict(mdl,X);
pred = (yfit > 0.5);
end

function pred = classf(X1train,X2train,X3train,ytrain,X1test,X2test,X3test)
Xtrain = table(X1train,X2train,X3train,ytrain);
Xtest = table(X1test,X2test,X3test, ...
    'VariableNames',{'Diastolic','Gender','Systolic'});
modelspec = 'Smoker ~ Diastolic + Gender + Systolic';
mdl = fitglm(Xtrain,modelspec,'Distribution','binomial');
yfit = predict(mdl,Xtest);
pred = (yfit > 0.5);
end
%}

%}

