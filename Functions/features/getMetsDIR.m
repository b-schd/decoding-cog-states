function [metric_matrices]= getMetsDIR(net, metrics)
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
  
    metrics={'SSIndex', 'sourceIndex', 'sinkIndex', 'strength',...
        'strengthIn', 'strengthOut', 'skewnessMet',...
        'kurtosisMet', 'clustering'};                                   % eigenstructure distribution %
end

[N, ~, T]= size(net);
config = reshape(net, N*N, T);


if ismember('skewnessMet', metrics), skewnessMet= skewness(config); end
if ismember('kurtosisMet', metrics), kurtosisMet= kurtosis(config); end


for t = 1:T
    
    A = net(:,:,t);
    
     if ismember('SSIndex', metrics)
        % Compute nodal ranks
        [~, i_col] = sort(sum(abs(A)));     col_rnk(i_col) = (1:N)./N;
        [~, i_row] = sort(sum(abs(A), 2));  row_rnk(i_row) = (1:N)./N;
         
        sinkIndex(:,t) = sqrt(2) - vecnorm([row_rnk; col_rnk] - [1;0])';
        sourceIndex(:,t) = sqrt(2) - vecnorm([row_rnk; col_rnk] - [0;1])';
        source_infl = A * sourceIndex(:,t);  
        sink_conn = A * sinkIndex(:,t);

        SSIndex(:,t) = sinkIndex(:,t) .* source_infl .* sink_conn; 
     end
     
     
     if ismember('strength', metrics)
       [inStr,outStr,strgth] = strengths_dir(A);
       strength(:,t) = strgth;
       strengthIn(:,t) = inStr;
       strengthOut(:,t) = outStr;
     end
     
     if ismember('clustering', metrics)
         clustering = clustering_coef_wd(abs(A));
     end
     
     
     
end

% Populate metric and state average structs
for m=metrics
    eval(sprintf('metric_matrices(1).%s=%s;',m{1},m{1}));
end
    
end