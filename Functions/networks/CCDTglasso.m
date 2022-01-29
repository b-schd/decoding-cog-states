function CCDTglasso(ddir, subjChLoc, subj, verbose, pdir, WstFilename)
%{
- ddir: data directory ('/Users/graceng/OneDrive/CCDTnewDat/CCDT/eeg/')
- subjChLoc: directory for file with information about channel location
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/'patient_loc.mat')
- subj: enter string name of the subject (eg, 'HUP069') or enter [] to
    process all subjects
- verbose: if true, prints a statement for each iteration of GLASSO (true or false)
- pdir (optional): processed data directory
    ('/Users/graceng/Documents/Med_School/Research/BCI/glasso_data/WstGLASSO/')
- WstFilename (optional): name of file storing Wst matrices
%}


% GLASSO parameters
numRho = 2; %100;
rhoRange = [0.01 0.1];
gamma = 0.5; %0.1 % goodness of fit index for EBIC; must be a value between [0, 1]
maxIt = 3; %1e2;


if nargin < 5
    pdir = [];
end
if nargin < 6
    WstFilename = 'WstGLASSO';
end

% check filepaths for saving data
if ~isempty(pdir)
    if ~isfolder(pdir)
        error('pdir does not exist.');
    end
end

% load processed data
[dataStruct, Nsubj] = preprocessCCDT(ddir, subjChLoc, subj);

% compute Wst using GLASSO to find effective connectivity
WstStruct = struct;
for isubj = 1:Nsubj
    datwin = dataStruct(isubj).datwin;
    Ntrl = dataStruct(isubj).Ntrl;
    Nch = dataStruct(isubj).Nch;
    Nsamp = dataStruct(isubj).Nsamp;
    subjName = dataStruct(isubj).subjName;
    
    % Compute an effective connectivity matrix for each trial
    allEffConn = zeros(Nch, Nch, Ntrl);
    convergences = zeros(Ntrl,1);
    rhoSet = linspace(rhoRange(1), rhoRange(2), numRho);
    for ii=1:Ntrl
        allTheta = zeros(numRho, Nch, Nch);
        allEBIC = zeros(numRho,1);
        allConverge = zeros(numRho,1);
        for jj = 1:numRho
            % Input to GraphicalLasso is a Nsamp x Nch matrix
            % Output is theta (dimensions are Nch x Nch)
            S = cov(transpose(squeeze(datwin(:,ii,:))));
            [Theta, ~, logL, converge] = graphicalLassoMatlab(S, rhoSet(jj), verbose, maxIt);
            allTheta(jj,:,:) = Theta;
            allConverge(jj) = converge;
            a_rho = nnz(Theta - diag(diag(Theta))); %a_rho is the number of non-zero off-diagonal elements in Theta(rho)
            allEBIC(jj) = -2*logL + a_rho*(log(Nsamp) + 4*gamma*log(Nch));
        end
        [~, idxMin] = min(allEBIC);
        bestTheta = squeeze(allTheta(idxMin,:,:));
        % standardize theta to generate an effective connectivity matrix
        allEffConn(:,:,ii) = -1. * bestTheta ./ sqrt(diag(bestTheta)) ./ sqrt(diag(bestTheta)).';
        convergences(ii) = allConverge(idxMin);
    end
    
    Wst = allEffConn;
    WstStruct(isubj).Wst = Wst;
    WstStruct(isubj).convergences = convergences;
    WstStruct(isubj).subjName = subjName;
    
    if ~isempty(pdir)
        save([pdir WstFilename, '.mat'],'WstStruct');
    end
end

end