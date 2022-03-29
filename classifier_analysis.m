% parameters
netType = 'AR_networks'; %name of directory with network metrics of interest
cleanTrials = true; %if true, remove all trials with negative reaction times
binarizeY = [33.3, 66.7]; %binarize trials according to these percentiles 
    %example: [40, 60] means that slowe trials are the lower 40% and 
    %fast trials are the highest 40%
netMetrics = {'strength', 'SSIndex', 'sourceIndex', 'sinkIndex'}; %network metric(s) of interest 
freqBands = {'ThetaAlpha', 'Beta', 'LowGamma','HighGamma'}; 
             %frequency band(s) of interest; use 'none' for no band passing
freqBandMap = containers.Map({'none', 'ThetaAlpha', 'Beta', 'LowGamma', ...
    'HighGamma'}, {[0 150], [3 12], [14 30], [35 55], [70 150]}); 
plrty = 1; %polarity can be 1 or 2 for channel location file
k_folds = 10; %number of folds for k-fold cross validation
nboot = 1000; %n for bootstrapping AUC values
metricNetType = 'Metrics_roi'; % Either "Metrics" or "MetricsRoi" 

% Load additional parameters
params = proj_config();

% Load SEEG location data 
load(params.subjChLoc);
db = CCDTdatabase;

for subj = {'HUP139', 'HUP154', 'HUP143', 'HUP165', 'HUP171', 'HUP179'}
    params.subj = subj{1}

subjInd = cellfun(@(x) strcmp(x,params.subj),db(:,1)); % find specified subj
if ~isempty(params.sess)
    subjInd2 = cellfun(@(x) strcmp(x,params.sess),db(:,2)); % find specified session
    subjInd = subjInd & subjInd2;
end
subjInd = find(subjInd);
if isempty(subjInd), error('subj/sess not found in database'); end
locInfo = patient_loc(plrty).session(subjInd);

% Load reaction time data
load(params.RTdata);
RT = RTstruct.(params.subj).RT;

% Option to clean (i.e., remove) trials where button was pressed before the cue
if cleanTrials
    trialIdxs = find(RT >= 0);
else
    trialIdxs = 1:length(RT);
end
Ntrl = length(trialIdxs);

% Set up structs to store results of analysis
stats = struct('preCue',struct,'preGo',struct);
for i = 1:length(freqBands)
    freqBand = freqBands{i};
    bands = freqBandMap(freqBand);
    stats.preCue.(freqBand) = {};
    stats.preGo.(freqBand) = {};

    % Load network metric data
    preCue = load(fullfile(params.ndir, params.subj, netType, ...
        sprintf('%s_net-ar_iev-1_bp-%d-%d.mat', params.subj, bands(1), bands(2))));
    preGo = load(fullfile(params.ndir, params.subj, netType, ...
        sprintf('%s_net-ar_iev-2_bp-%d-%d.mat', params.subj, bands(1), bands(2))));
    
    if ~ismember(metricNetType, fieldnames(preCue))
        disp(['preCue.', metricNetType, ' does not exist for this subject and band'])
        continue;  
    end
    
    for j = 1:length(netMetrics)
        netMetric = netMetrics{j};
        Nch = size(preCue.(metricNetType)(1).(netMetric),1);
        preCueMetric = ones(Ntrl,Nch,size(preCue.(metricNetType)(1).(netMetric),2));
        preGoMetric = ones(Ntrl,Nch,size(preGo.(metricNetType)(1).(netMetric),2));
        for k = 1:Ntrl
            preCueMetric(k,:,:) = preCue.(metricNetType)(trialIdxs(k)).(netMetric);
            preGoMetric(k,:,:) = preGo.(metricNetType)(trialIdxs(k)).(netMetric);
        end

        % Set up input data for classifier
        selectRT = RT(trialIdxs);
        RTprctile = prctile(selectRT, [binarizeY(1), binarizeY(2)]);
        fastIdxs = find(selectRT <= RTprctile(1));
        slowIdxs = find(selectRT > RTprctile(2));
        binarizeIdxs = union(slowIdxs, fastIdxs);
        tempY = -1.*ones(Ntrl,1);
        tempY(slowIdxs)=0;
        tempY(fastIdxs)=1;
        y = tempY(tempY~=-1);
        Xcue = preCueMetric(binarizeIdxs,:,end); 
        Xgo = preGoMetric(binarizeIdxs,:,end);
        %Note: "end" means that we only analyze network metrics from the last time bin 
        
        % Run classifier
        thisStats = runSVM(Xcue, Xgo, y, k_folds, nboot);
        stats.preCue.(freqBand).(netMetric).auc = thisStats{1}.auc;
        stats.preCue.(freqBand).(netMetric).scores = thisStats{1}.scores;
        stats.preCue.(freqBand).(netMetric).labels = thisStats{1}.labels;
        stats.preGo.(freqBand).(netMetric).auc = thisStats{2}.auc;
        stats.preGo.(freqBand).(netMetric).scores = thisStats{2}.scores;
        stats.preGo.(freqBand).(netMetric).labels = thisStats{2}.labels;
    end
end
if ~(params.odir == "")
    outFilepath = [params.odir, 'SVM_', params.subj, '_', date, '.mat'];
    save(outFilepath,'stats');
end
end

function stats = runSVM(Xcue, Xgo, y, k_folds, nboot)
    allX = {Xcue, Xgo};
    stats = {struct, struct}; %first struct stores stats about preCue data, 
                              %second struct stores stats about preGo data
    for m = 1:length(allX)
        thisX = allX{m};
        mdlSVM = fitclinear(thisX, y, 'Learner', 'svm', 'KFold', k_folds);
        [labels,scores] = kfoldPredict(mdlSVM);
        [~,~,~,auc] = perfcurve(y, scores(:,1), 1, 'NBoot', ...
            nboot, 'boottype', 'cper');
        if auc(1) < 0.5 %if AUC < 0.5, flip labels for positive and negative classes
            [~,~,~,auc] = perfcurve(y, scores(:,1), 0, 'NBoot', ...
                nboot, 'boottype', 'cper');
        end
        stats{m}.auc = auc;
        stats{m}.scores = scores;
        stats{m}.labels = labels;
    end
end

%%%%%%%%%%%%%%%%%%%%% OLD CODE -- needs cleaning %%%%%%%%%%%%%%%%%%%%%%%

% Multiple linear regression
%{
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
%}

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

