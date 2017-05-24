function rdmSet = rdm_computeMahalRDMset(betas,resids)
%% RDM_COMPUTEMAHALRDM
% 
% computes mahalanobis distance RDMset 
%
% betas  = electrode-time-condition matrix
% resids = electrode-time-trial matrix
%
% (c) Timo Flesch, 2016

rdmSet = [];
% we want a time-condition-condition rdm-set:
rdmSet = zeros(size(betas,2),size(betas,3),size(betas,3));

% iterate through all time points
    for timePoint = 1:(size(betas,2))    
        respMat = squeeze(mean(betas(:,timePoint,:),2));
        residMat = squeeze(mean(resids(:,timePoint,:),2));
        rdmSet(timePoint,:,:) = squareform(pdist(respMat','mahalanobis',covdiag(residMat')));
    end
   
end