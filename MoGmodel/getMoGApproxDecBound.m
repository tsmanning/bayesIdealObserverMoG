function [decisionMask] = getMoGApproxDecBound(support1,support2,priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2)

% Given a 2D likelihood function, what's the probability of responding s2
% instad of s1? Approximated numerically with Eqn. 18 & 36

supportY = fliplr(support2);

% Get number of support elements
suppEl = numel(support1);

decisionMask     = nan(suppEl,suppEl);

% Calculate whether each measurement pair would be above or below decision
% boundary
for ii = 1:suppEl
    for jj = 1:suppEl
        
        thisS1measurement = support1(ii);
        thisS2measurement = supportY(jj);
        
        % Collect modified weights and means and the shrinkage factors that define
        % posteriors for each measurement pair
        [alphas1,muTildes1,w1] = getMogPostPars(priorWs,priorSigs,priorMus,sig1,thisS1measurement);
        [alphas2,muTildes2,w2] = getMogPostPars(priorWs,priorSigs,priorMus,sig2,thisS2measurement);
        
        [~,~,wTildes1EV] = getMogPostPars(priorWs,priorSigs,priorMus,sig1,mu1);
        [~,~,wTildes2EV] = getMogPostPars(priorWs,priorSigs,priorMus,sig2,mu2);
        
        post1 = wTildes1EV'*(alphas1*thisS1measurement + muTildes1);
        post2 = wTildes2EV'*(alphas2*thisS2measurement + muTildes2);
        
        decisionMask(jj,ii) = post2 > post1;
    end
end


end