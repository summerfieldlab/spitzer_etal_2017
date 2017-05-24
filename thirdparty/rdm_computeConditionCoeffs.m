    function results = rdm_computeConditionCoeffs(eegMat,behavVect)
    %% RDM_COMPUTECONDITIONCOEFFS
    %
    % regresses the eeg signal on the trial-by-condition design matrix for all electrodes and points in time
    %
    % eegMat = electrode-by-time-by-trial
    % behavVect = trial-by-1 (indicating the condition on each trial)
    %
    % (c) Timo Flesch, 2016
     
    % set up design matrix
    conditions = unique(behavVect)';
    dmat = zeros(size(behavVect,1),length(conditions));
    for c = 1:length(conditions)
        dmat(behavVect==conditions(c),c) = 1;   
    end
    % add constant term
    dmat = [ones(size(dmat,1),1) dmat];
    
    results = struct(); 
    results.betas   = zeros(size(eegMat,1),size(eegMat,2),length(conditions));
    results.resids  = zeros(size(eegMat));
     
    % for all electrodes and time points: run regression
    for el = 1:size(eegMat,1)
        %el
        for t = 1:size(eegMat,2)
           [betas,~,resid] = regress(zscore(squeeze(eegMat(el,t,:))),dmat);
           results.betas(el,t,:) = betas(2:end); % drop constant term
           results.resids(el,t,:) = resid;
        end       
    end
    end

