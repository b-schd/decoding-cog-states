function [betahat,beta,Rsquared,RMSE] = bootstrapreg(x,y,N,alpha)
%% Conducting bootstrap linear regression model of y = b0 + b1*x1 + b2*x2 + ... + bn*xn + error
% INPUT
% x: explanatory variable(s)
% y: response variable
% N: number of bootstrping simulations (e.g. 1000)
% alpha: significant level (e.g., 0,05, 0.1)
%
% OUTPUT
% betahat: a structure containing the mean and confidence interval of model
% parameters
% beta: random samples of model parameters generated from bootstrap simulation
% Rsquared: random samples of Rsquared generated from bootstrap simulation
% RMSEP: random samples of RMSEP generated from bootstrap simulation
%% Copyright (C): May 15, 2020 
% Created by Shiyong Yu
% School of Geographical Sciences
% Jiangsu Normal University
% China 
% E-mail: syu@jsnu.edu.cn 
% To cite this file, this would be an appropriate format:
% Yu, S.-Y., (2020). Bootstrapreg: A MATLAB program for bootstrapping a
% multivariate linear regression model. 
% (https://www.mathworks.com/matlabcentral/fileexchange/), 
% MATLAB Central File Exchange. Retrieved May 18, 2020.
%% Boostrap simulation 
y = y(:);
if size(x,2) == length(y)
   x = x';
end
M = size(x,2) + 1; %number of model parameters
beta = zeros(N,M);
betapv = zeros(N,M);
RMSE = zeros(N,1);
Rsquared = zeros(N,1);
[yhat,t] = tresidual(x,y);              %initialization
for i = 1:N
    t_star = t(randperm(length(t)));    %resampling the residual with replacement
    y = yhat + t_star;                  %updating the response variable 
    lm = fitlm(x,y,'linear');           %fitting the linear model 
    est = lm.Coefficients.Estimate;     %extracting the estimate of model parameters
    p = lm.Coefficients.pValue;         %extracting the P value of model parameter
    RMSE(i) = lm.RMSE;                  %extracting the RMSE
    Rsquared(i) = lm.Rsquared.ordinary; %extracting R2
    beta(i,:) = est';
    betapv(i,:) = p';
end
%% Calcualting mean and confidence interval of model parameters
betahat.mean = mean(beta,1)';
betahat.CI = zeros(M,2);
for j = 1:M
    betahat.CI(j,:) = prctile(beta(:,j),abs([0,100]-(100-100*(1-alpha/2))));
end
%% Calculating p values of model parameters
betahat.pvalue = mean(betapv,1)';
%% Plotting results
figure(1)
for j = 1:M
    subplot(1,M,j); histogram(beta(:,j),50);
    hold on;
    ylim = get(gca,'YLim');
    plot(betahat.CI(j,1)*[1,1],ylim,'r-','LineWidth',1); %adding the lower 100*(1-alpha) confidence level
    plot(betahat.CI(j,2)*[1,1],ylim,'r-','LineWidth',1); %adding the upper 100*(1-alpha) confidence level
    xlabel(['{\beta}_' num2str(j-1)]); 
    ylabel('Number of bootstrap samples');
end
figure(2)
subplot(1,2,1); histogram(Rsquared,50);
xlabel('R^{2}'); 
ylabel('Number of bootstrap samples');
subplot(1,2,2); histogram(RMSE,50);
xlabel('RMSE'); 
ylabel('Number of bootstrap samples');
%% Saving results
results = [beta Rsquared RMSE];
header = cell(1,M+2);
for j = 1:M
    header{j} = strcat('beta_', num2str(j-1));
end    
header{M+1} ='Rsquared';
header{M+2} = 'RMSE';
fmt = repmat('%12g ',1,M+2); 
fmt = [fmt '\n'];
fid = fopen('bootstrap_results','wt');
fprintf(fid,'%12s ',header{1:end});
fprintf(fid,'\n');
for i = 1:N
    fprintf(fid,fmt,results(i,:));
end    
fclose(fid);
return;
%%
function [yhat,t] = tresidual(x,y)
%Calculating the studentized residuals
%INPUT
%x: a vector containing the explanatory variable(s)
%y: an array containing the response variable
%OUTPUT
%t: an array containing the studentized residuals
%yhat: an array containing the predicted response variable
%%
y = y(:);
if size(x,2) == length(y)
    x = x';
end
%Building up the X matrix
[M,N] = size(x);
X = ones(M,N+1);
X(:,2:end) = x;
%Bulding up the hat matrix
H = X*((X'*X)\X');
h = diag(H); 
%Generating predictions
lm = fitlm(x,y,'linear');
yhat = predict(lm,x);
%Calculating the studentized residual
epsilon = y - yhat;
sigma2 = (epsilon'*epsilon)/(M-N-1);
sigma = sqrt(sigma2);
t = epsilon./(sigma.*sqrt(1-h));
return;