%% Define load and save directories
loadDir = '/Users/graceng/Documents/Med_School/Research/BCI/CCDT mat data';
saveDir = '/Users/graceng/Documents/Med_School/Research/BCI/CCDT mat data';
dataDir = '/allSubjBehavAndNetworkData';

%% Define analysis parameters
periods = {'precueone' 'precue' 'DTone' 'DTtwo' 'DTthree'};
minRT = 0; %the minimum reaction time in order for a trial to qualify as 
            % "no-error"
DTThresh = 50; %trials with delay times within this threshold are considered 
                %in this analysis
                
RToTdata = 'r_RToT_allNoE_sI0';
netStrType = "norm"; %"str";
%typeTrial = "E"; %if typeTrial is "E", then we analyze the error trials
typeTrial = "noE"; %if typeTrial is "noE", the we analyze the non-error trials
sITypes = [0]; %[0,1,2]

numTimeBins = 4; %average HG Net Str according to these time bins before analysis
avgTimeBins = false;
avgScoreBins = 10;

%DT = 1500; %the delay time of the trials to be considered in this analysis
%typeDT = "min"; %if typeDT is 'min', then we consider trials with DTs greater 
                    %than or equal to DT-DTThresh
%DT = 500;
%typeDT = "max"; %if typeDT is 'max', then we consider trials with DTs less
                    %than or equal to DT+DTThresh
DT = [1500, 500];
typeDT = ["min", "max"];                

freqInd = 3; %the index for the frequency band data 
                %(1 for theta, 2 for high gamma, 3 for low gamma)
freqType = 'LG';
%freqInd = 1;
%freqType = 'theta';
%freqInd = 2;
%freqType = 'HG'


%% run analysis
load([loadDir dataDir], 'allSubjBehaviorStruct', 'allSubjNetMetrics',...
    'sortedLearners','subjOrder',RToTdata);
RToT = eval(RToTdata);
ind = [];
minTrials = 0;
qualSubjects = []; 
if length(DT) > 1 
    for i=1:length(DT)
        [thisInd, thisMinTrials, thisQualSubjects] = ...
            filterInputData(allSubjBehaviorStruct, minRT, DTThresh, ...
            typeTrial, sITypes, DT(i), typeDT(i), subjOrder);
        if isempty(ind)
            ind = thisInd;
        else
            for j=1:length(thisInd)
                ind{j} = ind{j}|thisInd{j};
            end
        end
        minTrials = minTrials + thisMinTrials;
        qualSubjects = union(qualSubjects, thisQualSubjects);
    end
else
    [ind, minTrials, qualSubjects] = filterInputData(allSubjBehaviorStruct, ...
        minRT, DTThresh, typeTrial, sITypes, DT, typeDT, subjOrder);
end
if ~avgTimeBins
    [X, y] = setInputsOutputs(allSubjNetMetrics, RToT, periods, netStrType, ...
        freqInd, ind, qualSubjects, minTrials); 
else
    [X, y] = setAvgInputsOutputs(allSubjNetMetrics, RToT, periods, ...
        netStrType, freqInd, ind, qualSubjects, numTimeBins);
    minTrials = numTimeBins;
end
[n, p] = size(X);
if n/p < 8
    %Use the generated graph to select the number of components for PLS and PCR
    plotCrossVal(X,y);
    numPLS = 3;
    numPCR = 7;
    orderLearnQualSubj = intersect(sortedLearners, qualSubjects, 'stable');
    [sortedSubj, sortIdxs] = sort(orderLearnQualSubj);
    orderLearnIdxs(sortIdxs) = 1:length(qualSubjects);
    PLSanalysis(X, y, orderLearnIdxs, qualSubjects, numPLS, freqType, DT, ...
        typeTrial, sITypes, periods, minTrials, avgScoreBins, subjOrder);
    %PCAanalysis(X,y,numPCR,periods,minTrials);
else
    linReg(X, y, DT, typeTrial, sITypes, periods, minTrials);
end

%% functions to run the analysis
function [ind, minTrials, qualSubjects] = filterInputData(...
    allSubjBehaviorStruct, minRT, DTThresh, typeTrial, sITypes, DT, ...
    typeDT, subjOrder)
    
    if typeDT == "min"
        indNoE = cellfun(@(x) x.vRT>minRT & x.vDT>(DT-DTThresh), ...
            allSubjBehaviorStruct, 'UniformOutput',false);
    elseif typeDT == "max"
        indNoE = cellfun(@(x) x.vRT>minRT & x.vDT<(DT+DTThresh), ...
            allSubjBehaviorStruct, 'UniformOutput',false);
    else
        error('typeDT must either be "min" or "max".');
    end

    if typeTrial == "E"
        ind = cellfun(@(x) ~x, indNoE);
    elseif typeTrial == "noE"
        ind = indNoE;
    else
        error('typeTrial must either be "E" or "noE".');
    end
    
    qualSubjects = [];
    minTrials = NaN;
    for i=1:length(subjOrder)
        if isfield(allSubjBehaviorStruct{i},"sI")
            sIidxs = zeros(length(ind{i}), 1);
            for j=1:length(sITypes)
                sIidxs = sIidxs|allSubjBehaviorStruct{i}.sI==sITypes(j);
            end
            ind{i} = ind{i}&sIidxs;
        else
            if ~ismember(0,sITypes)
                ind{i} = NaN;
            end
        end
        if ~isnan(ind{i})
            if any(ind{i}) 
                qualSubjects = [qualSubjects, i];  %#ok<AGROW>
                if isnan(minTrials)
                    minTrials = sum(ind{i});
                else
                    if minTrials > sum(ind{i})
                        minTrials = sum(ind{i});
                    end
                end
            end
        end
    end
end

function [X, y] = setInputsOutputs(allSubjNetMetrics, RToT, periods, ...
    netStrType, freqInd, ind, qualSubjects, minTrials)
    netStr = zeros(length(qualSubjects), minTrials, length(periods));
    transformedRToT = zeros(length(qualSubjects), 1);
    for i=1:length(qualSubjects)
        idx = qualSubjects(i);
        sampleInds = zeros(sum(ind{idx}),1);
        %sampleInds(randsample(sum(ind{idx}), minTrials)) = 1; 
        %Use this line instead of the next if you prefer to select trials 
        %randomly instead of at regularly-spaced intervals
        sampleInds(round(linspace(1,sum(ind{idx}),minTrials))) = 1;
        sampleInds = logical(sampleInds);
        for j=1:length(periods)
            if netStrType == "str"
                dataLabel = ['str_' periods{j}];
            elseif netStrType == "norm"
                dataLabel = ['normStr_' periods{j}];
            else
                error('netStrType must either be "str" or "norm".');
            end
            sampleMean = mean(allSubjNetMetrics{idx}.(dataLabel)(...
                (ind{idx}),:,freqInd),2);
            netStr(i,:,j) = sampleMean(sampleInds);
        end
        transformedRToT(i) = 1. - RToT(idx);
    end
    X = reshape(netStr, length(qualSubjects), minTrials*length(periods));
    [n,~] = size(X);
    y = reshape(transformedRToT,n,1);
end

function [X, y] = setAvgInputsOutputs(allSubjNetMetrics, RToT, periods, ...
    netStrType, freqInd, ind, qualSubjects, numTimeBins)
    netStr = zeros(length(qualSubjects), numTimeBins, length(periods));
    transformedRToT = zeros(length(qualSubjects), 1);
    for i=1:length(qualSubjects)
        idx = qualSubjects(i);
        for j=1:length(periods)
            if netStrType == "str"
                dataLabel = ['str_' periods{j}];
            elseif netStrType == "norm"
                dataLabel = ['normStr_' periods{j}];
            else
                error('netStrType must either be "str" or "norm".');
            end
            sampleMean = mean(allSubjNetMetrics{idx}.(dataLabel)(...
                (ind{idx}),:,freqInd),2);
            r = diff(fix(linspace(0, sum(ind{idx}), numTimeBins+1)));
            binnedSampleMean = mat2cell(sampleMean, r, 1);
            for k=1:numTimeBins
                netStr(i,k,j) = mean(binnedSampleMean{k});
            end
        end
        transformedRToT(i) = 1. - RToT(idx);
    end
    X = reshape(netStr, length(qualSubjects), numTimeBins*length(periods));
    [n,~] = size(X);
    y = reshape(transformedRToT,n,1);
end

function plotCrossVal(X,y)
    [n,~] = size(X);
    if n > 100
        numCrossVal = 10;
    elseif n > 10
        numCrossVal = 5;
    else
        return
    end
    [~,~,~,~,~,~,PLSmsep] = plsregress(X,y,10,'CV',numCrossVal);
    PCRmsep = sum(crossval(@pcrsse,X,y,'KFold',10),1) / n;
    figure
    plot(0:10,PLSmsep(2,:),'b-o',0:10,PCRmsep,'r-^');
    xlabel('Number of components');
    ylabel('Estimated Mean Squared Prediction Error');
    legend({'PLSR' 'PCR'},'location','NE');
    title('Selecting the Number of Components');
end

function PLSanalysis(X, y, orderLearnIdxs, qualSubjects, numPLS, freqType, ... 
    DT, typeTrial, sITypes, periods, minTrials, avgScoreBins, subjOrder)
    [Xloadings,Yloadings,Xscores,~,betaPLS,pctVar,mse,stats] = ...
        plsregress(X,y,numPLS);
    [n,p] = size(X);
    yfitPLS = [ones(n,1) X]*betaPLS;
    % Calculate normalized PLS weights
    W0 = bsxfun(@rdivide,stats.W,sqrt(sum(stats.W.^2,1)));
    % Calculate the product of summed squares of XS and YL
    sumSq = sum(Xscores.^2,1).*sum(Yloadings.^2,1);
    % Calculate VIP scores for NCOMP components
    vipScores = sqrt(size(Xloadings,1) * sum(bsxfun(@times,sumSq,W0.^2),2)...
        ./ sum(sumSq,2));
    plotVIPScores(periods, vipScores, minTrials, freqType, DT, typeTrial, ...
        sITypes, avgScoreBins);
    plotPLSError(y, yfitPLS, orderLearnIdxs, qualSubjects, numPLS, DT, ...
        typeTrial, sITypes, mse, pctVar, subjOrder);
    plotPLSComponents(p, stats, periods, minTrials, numPLS);

    TSS = sum((y-mean(y)).^2);
    RSS_PLS = sum((y-yfitPLS).^2);
    rsquaredPLS = 1 - RSS_PLS/TSS;
    fprintf('rsquared PLS is %f\n', rsquaredPLS);
end

function plotVIPScores(periods, vipScores, minTrials, freqType, DT, ...
    typeTrial, sITypes, avgScoreBins)
    figure
    hold on
    for i=1:length(periods)
        plot((i-1)*minTrials+1:i*minTrials, ...
            vipScores((i-1)*minTrials+1:i*minTrials),'.-',... 
            'MarkerSize',16)
    end
    xlabel('Variable');
    xticks((0:length(periods)-1)*minTrials+1);
    xticklabels(periods);
    ylabel('PLS Variable in Projection Score');
    title(join(['Distribution of PLS VIP Scores for' freqType 'Net Str' ...
        num2str(DT,'%d') 'ms' typeTrial 'Trials: Session(s)' sITypes]));
    yline(1.,'--r');
    dim = [.77 .82 .1 .1];
    str = '*Scores greater than 1 are relevant';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
    ax = gca;
    ax.XGrid = 'on';
    
    if avgScoreBins > 0
        figure
        hold on
        r1 = diff(fix(linspace(0, length(vipScores), length(periods)+1)));
        scoresByPeriod = mat2cell(vipScores, r1, 1);
        for i=1:length(periods)
            r2 = diff(fix(linspace(0, length(scoresByPeriod{i}), ...
                avgScoreBins+1)));
            binnedScores = mat2cell(scoresByPeriod{i}, r2, 1);
            means = cellfun(@(x) mean(x), binnedScores);
            stdDevs = cellfun(@(x) std(x), binnedScores);
            errorbar((i-1)*avgScoreBins+1:i*avgScoreBins, ...
            means, stdDevs, '.-', 'MarkerSize',16);
        end
        xlabel('Variable');
        xticks((0:length(periods)-1)*avgScoreBins+1);
        xticklabels(periods);
        ylabel('Averaged PLS Variable in Projection Scores');
        title(join(['Distribution of Averaged PLS VIP Scores for' ...
            freqType 'Net Str' num2str(DT,'%d') 'ms ' typeTrial ...
            'Trials: Session(s)' sITypes]));
        yline(1.,'--r');
        dim = [.77 .82 .1 .1];
        str = '*Scores greater than 1 are relevant';
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
        ax = gca;
        ax.XGrid = 'on';
    end
end

function plotPLSError(y, yfitPLS, orderLearnIdxs, qualSubjects, numPLS, ...
    DT, typeTrial, sITypes, mse, pctVar, subjOrder)

    figure
    residuals = y - yfitPLS;
    sortedResiduals = residuals(orderLearnIdxs);
    stem(sortedResiduals)
    xlabel('Observations (sorted by learning rate)');
    ylabel('Residual');
    xticks(1:length(y));
    subjIdxs = qualSubjects(orderLearnIdxs);
    xticklabels(subjOrder(subjIdxs));
    xtickangle(45)
    title(join(['Residuals of PLS Model for ' num2str(DT,'%d') 'ms ' ...
        typeTrial ' Trials: Session(s) ' sITypes]));

    figure
    plot(y,yfitPLS,'bo');
    xlabel('Observed Response');
    ylabel('Fitted Response');
    title("Comparing Observed and Fitted Responses for PLS with " + ...
        num2str(numPLS) + " Components")

    figure
    yyaxis left
    plot(0:numPLS,mse(1,:));
    ylabel ('MSE Predictors');
    yyaxis right
    plot(0:numPLS,mse(2,:));
    ylabel('MSE Response');
    xticks(1:numPLS);
    xlabel('Number of Components')
    title('Mean-Squared Error of the PLS Model')


    figure
    yyaxis left
    plot(1:numPLS,pctVar(1,:));
    ylabel ('% of Predictor Variance Explained');
    yyaxis right
    plot(1:numPLS,pctVar(2,:));
    ylabel('% of Response Variance Explained');
    xticks(1:numPLS);
    xlabel('Number of Components')
    title('Percentage Variance Explained by the PLS Model')
end

function plotPLSComponents(p, stats, periods, minTrials, numPLS)
    figure
    plot(1:p,stats.W,'-');
    xlabel('Variable');
    xticks((0:length(periods)-1)*minTrials+1);
    xticklabels(periods);
    ylabel('PLS Weight');
    legendStrings = cellfun(@(c) "Component " + num2str(c), num2cell(1:numPLS));
    legend(legendStrings, 'location','NW');
    title('Weight Distribution of PLS Components');

    for i=1:numPLS
        figure
        plot(1:p,stats.W(:,i),'-');
        xlabel('Variable');
        xticks((0:length(periods)-1)*minTrials+1);
        xticklabels(periods);
        ylabel('PLS Weight');
        legend({"Component " + num2str(i)});
        title("Weight Distribution of PLS Component " + num2str(i));
    end
end

function PCAanalysis(X, y, numPCR, periods, minTrials)
    [n,p] = size(X);
    [PCALoadings,PCAScores,~] = pca(X,'Economy',false);
    betaPCR = regress(y-mean(y), PCAScores(:,1:numPCR));
    betaPCR = PCALoadings(:,1:numPCR)*betaPCR;
    betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];
    yfitPCR = [ones(n,1) X]*betaPCR;
    
    figure
    plot(y,yfitPCR,'r^');
    xlabel('Observed Response');
    ylabel('Fitted Response');
    title("Comparing Observed and Fitted Responses for PCR with " + ...
        num2str(numPCR) + " Components");
    
    TSS = sum((y-mean(y)).^2);
    RSS_PCR = sum((y-yfitPCR).^2);
    rsquaredPCR = 1 - RSS_PCR/TSS;
    fprintf('rsquared PCR is %f\n', rsquaredPCR);    
    
    figure
    plot(1:p,PCALoadings(:,1:numPCR),'-');
    xlabel('Variable');
    xticks((0:length(periods)-1)*minTrials+1);
    xticklabels(periods);
    ylabel('PCA Loading');
    legendStrings = cellfun(@(c) "Component " + num2str(c), num2cell(1:numPCR));
    legend(legendStrings, 'location','NW');
    title('Weight Distribution of PCA Components');

    for i=1:numPCR
        figure
        plot(1:p,PCALoadings(:,i),'-');
        xlabel('Variable');
        xticks((0:length(periods)-1)*minTrials+1);
        xticklabels(periods);
        ylabel('PCA Loading');
        legend({"Component " + num2str(i)});
        title("Weight Distribution of PCA Component " + num2str(i));
    end
end

function linReg(X, y, DT, typeTrial, sITypes, periods, minTrials)
    mdl = fitlm(X,y);
    beta = mdl.Coefficients.Estimate;
    
    figure
    hold on
    for i=1:length(periods)
        plot((i-1)*minTrials+1:i*minTrials, ...
            beta((i-1)*minTrials+1:i*minTrials),'.-',... 
            'MarkerSize',16)
    end
    xlabel('Variable');
    xticks((0:length(periods)-1)*minTrials+1);
    xticklabels(periods);
    ylabel('Linear Regression Weights');
    title(join(['Distribution of Linear Regression Weights for ' ...
        num2str(DT,'%d') 'ms ' typeTrial ' Trials: Session(s) ' sITypes]));
    yline(1.,'--r');
end