function [pFxn,decisionMask] = calcMoGPFxn_Numeric(support1,support2,priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2,plotOn)

% Given a 2D likelihood function, what's the probability of responding s2
% instad of s1? Approximated numerically with Eqn. 18 & 36

supportY = fliplr(support2);

% Get number of support elements
suppEl = numel(support1);

% Generate 2D likelihood function
like1 = normpdf(support1,mu1,sig1);
like2 = normpdf(supportY,mu2,sig2);

like2D = like2'*like1;
like2D = like2D/sum(like2D(:));

decisionMask     = nan(suppEl,suppEl);
decisionBoundary = decisionMask;

% Calculate whether each measurement pair would be above or below decision
% boundary
for ii = 1:suppEl
    for jj = 1:suppEl
        
        thisS1measurement = support1(ii);
        thisS2measurement = supportY(jj);
        
        % Collect modified weights and means and the shrinkage factors that define
        % posteriors for each measurement pair
        [alphas1,muTildes1,wTildes1] = getMogPostPars(priorWs,priorSigs,priorMus,sig1,thisS1measurement);
        [alphas2,muTildes2,wTildes2] = getMogPostPars(priorWs,priorSigs,priorMus,sig2,thisS2measurement);
        
        post1 = wTildes1'*(alphas1*thisS1measurement + muTildes1);
        post2 = wTildes2'*(alphas2*thisS2measurement + muTildes2);
        
        decisionMask(jj,ii)     = post2 > post1;
        decisionBoundary(jj,ii) = post2 == post1;
    end
end

% Mask likelihood with decision matrix and integrate above boundary
pFxn = sum(decisionMask(:).*like2D(:)) + 0.5*sum(decisionBoundary(:).*like2D(:));

if plotOn
   
    f1 = figure;
    f1.Position = [100 100 650 650];
    hold on;
    
    imagesc(support1,supportY,decisionMask.*like2D);
    set(gca,'plotboxaspectratio',[1 1 1],'ydir','reverse','xlim',[support1(1) support1(end)],...
            'ylim',[support2(1) support2(end)],'fontsize',20,'xtick',[],'ytick',[]);
    xlabel('m_{1}');
    ylabel('m_{2}');
    
end


end