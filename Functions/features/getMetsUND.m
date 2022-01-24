function [metric_matrices]=getMetsUND(Network, metrics, thresh)
% Get metrics for UNDIRECTED (symmetric) adjacency matrices
% INPUTS:
%   Network: A struct with fields pcm and config_pcm holding partial corr. matrices
%   metrics: a string array with the metric names to calculate:
%       'globalCtrl', 'aveCtrl', 'modalCtrl', 'pmodalCtrl', 'tModalCtrl', 
%       'strength', 'strenghtPos', 'strenghtNeg', 'sekwnessMet',
%       'kurtosisMet', degree', 'clustering3'
%
% OUTPUTS:
%   metric_matrices: Struct with fields corresponding to the calculated
%                    metrics. 


metric_matrices= struct(); 

if nargin < 2
    thresh=0.15;     % Threshold for transient/persistent mode selection
    metrics={'globalCtrl', 'aveCtrl', 'modalCtrl', 'pModalCtrl', 'tModalCtrl',...
       'strength', 'strengthPos', 'strengthNeg', ...        % Network metrics %
       'skewnessMet', 'kurtosisMet', 'degree', 'clustering3'}; %, ... % connection distribution %
       %'spreadEig', 'skeweig'};                             % eigenstructure distribution %
end
dt=1;  

% only compute for beta= 0.01


% match ID

Net=Network(1); 

% Initialize metric structs
% metric_matrices(i_set).ID= p.ID; State_metrics(i_set).ID=p.ID;    
% metric_matrices(i_set).type= p.type; State_metrics(i_set).type=p.type; 
% metric_matrices(i_set).block=p.block; State_metrics(i_set).block=p.block;
% metric_matrices(i_set).beta=p.beta; State_metrics(i_set).beta=p.beta;

config=  Net.config_pcm;
pcm=     Net.pcm;  
[N, ~, T]=       size(pcm);

% Initialize average, modal metrics
[globalCtrl, aveCtrl, modalCtrl, pModalCtrl, tModalCtrl,...
   strength, strengthPos, strengthNeg, ...
   degree, clustering3]= deal(zeros(N, T));

if ismember('skewnessMet', metrics), skewnessMet= skewness(config); end
if ismember('kurtosisMet', metrics), kurtosisMet=kurtosis(config); end

for t=1:T
   if mod(t,10)==0; fprintf('%s', '.'); end

   %%%%%%%%%%%%%%%%%%%%%%%%
   %%% Control  Metrics %%%
   %%%%%%%%%%%%%%%%%%%%%%%%

   % Normalize A
    A= pcm(:,:,t);
    normA= A./(1+svds(A,1))-eye(N);

    % Global Control
   if ismember('globalCtrl', metrics)
       for n=1:N
           B=eye(N).*1e-6; B(n,n)=1;                % set control points
           sys= ss(normA, B, [], []);
           gramian= gram(sys, 'c');
           globalCtrl(n,t)=min(eig(gramian)); 
       end
   end

   % Avg Control
   if ismember('aveCtrl', metrics)
       sys=ss(normA,eye(N), [], []); 
       gramian=gram(sys, 'c');
       aveCtrl(:,t)= diag(gramian);
   end

   %Modal Control--- Note: norm happens in function
   if ismember('modalCtrl', metrics), modalCtrl(:,t)= modal_control(A, dt); end 

   % pModal and tModal 
   % Note: norm happens in func, find best threshold using "assessModalThresh.m"
  if ismember('pModalCtrl', metrics)
    [pModalCtrl(:,t), tModalCtrl(:,t)]= modal_persist_trans(A, thresh, dt);
  end


   %%%%%%%%%%%%%%%%%%%%%%%%
   %%% Network Metrics %%%
   %%%%%%%%%%%%%%%%%%%%%%%%
   if ismember('degree', metrics), degree(:,t)=degrees_und(pcm(:,:,t)); end

   % Calculating strength as sum of pos and abs(neg)
   if ismember('strength', metrics)
       [nPos,nNeg,vpos,vneg]=strengths_und_sign(pcm(:,:,t));
       strength(:,t)=(nPos'-nNeg');
       strengthPos(:,t)=nPos;
       strengthNeg(:,t)=nNeg;
   end

   % Calculating signed clustering using Perugini's signed extension
   if ismember('clustering3', metrics)
       [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(pcm(:,:,t),3);
       clustering3(:,t)=C_pos; 
   end

   %%%%%%%%%%%%%%%%%%%%%%%%
   %%%  Eigenstructure  %%%
   %%%%%%%%%%%%%%%%%%%%%%%%
   %spreadEig(:,t)=eig(normA);

end

% Populate metric and state average structs
for m=metrics
    eval(sprintf('metric_matrices(1).%s=%s;',m{1},m{1}));
end
    
end



% Additional strength metrics:
% for i_set=1:length(Metric_matrices)
%     %Metric_matrices(i_set).strengthAbs=abs(Metric_matrices(i_set).strengthPos)+ abs(Metric_matrices(i_set).strengthNeg)
%     %Metric_matrices(i_set).PosRatio=sum(Metric_matrices(i_set).strengthPos)./sum(Metric_matrices(i_set).strengthAbs)
% end

