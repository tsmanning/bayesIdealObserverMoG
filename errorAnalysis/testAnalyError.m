% Testing script to demonstrate error in analytical approximation vs
% numerical one

clear all; 
close all;

colorMat = colororder;

%% Setup simulation pars

% Define number of prior/likelihood combinations
numSims = 5000;
% numSims = 100;

% Define bounds for stimuli/measurements
stimBnds    = [-1 1];

% Define bounds for observer priors
% priorConsts = 'gtLike';
% priorConsts = 'longTail';
priorConsts = 'bimodal';
% priorConsts = 'bimodalLikeOut';
% priorConsts = 'none';

% Set prior scale factor (e.g. 1 - prior components are 1-2x wider than
% likelihoods, 2 - prior components are 1-4x wider, etc)
prScale     = 3;

% Set max number of components in prior
% maxComps    = 5;
maxComps    = 2;


%% Run simulations

simMat = struct;

for ii = 1:numSims
    
    if mod(ii,10) == 0
        disp(['Running simulation: ',num2str(ii),'/',num2str(numSims)]);
    end
    
    % Randomly select stimuli, define random measurement distributions
    x1   = rand*stimBnds(2)*2 - stimBnds(2);
    mu1  = x1;
    sig1 = rand*stimBnds(2)*0.5;
    
    x2   = rand*stimBnds(2)*2 - stimBnds(2);
    mu2  = x2;
    sig2 = rand*stimBnds(2)*0.5;
    
    % Define a random prior    
    switch priorConsts
        case 'gtLike'
            % Constrain prior s.t. it's narrowest component is wider than
            % the widest likelihood
            numComps  = randi(maxComps-1) + 1;
    
            priorWs   = rand(1,numComps);
            priorWs   = priorWs/sum(priorWs);
            priorMus  = rand(1,numComps)*stimBnds(2) - stimBnds(2)/2;
            priorSigs = max([sig1 sig2])*(1 + rand(1,numComps)*prScale);
            
        case 'longTail'
            % Constrain prior s.t. it has a long tail and means are both
            % zero
            numComps  = randi(maxComps-1) + 1;
    
            priorWs   = rand(1,numComps);
            priorWs   = priorWs/sum(priorWs);
%             priorMus  = [0 rand(1,numComps-1)*stimBnds(2) - stimBnds(2)/2];            
            priorMus  = zeros(1,numComps);
%             priorSigs = [1.1*max([sig1 sig2]) max([sig1 sig2])*(1.25 + rand(1,numComps-1)*prScale)];
            priorSigs = [1.1*max([sig1 sig2]) max([sig1 sig2])*(1.75 + rand(1,numComps-1)*prScale)];
            
        case 'bimodal'
            % Constrain prior to be bimodal
            numComps  = randi(maxComps-1) + 1;
    
            priorWs   = 0.3*ones(1,numComps)+rand(1,numComps);
%             priorWs   = ones(1,numComps);
            priorWs   = priorWs/sum(priorWs);
            
            priorSigs = max([sig1 sig2])*(1 + rand(1,numComps)*prScale*0.125);
            
            if numComps < 4
%                 priorMus  = [rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
%                              rand(1,numComps-1)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
                         
                priorMus(1)  = rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5;
                priorMus(2)  = priorMus(1) + 3*max(priorSigs) + rand(1)*stimBnds(2)*0.4;
            else
                priorMus  = [rand(1,2)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
                             rand(1,numComps-2)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
            end
            
            
        case 'bimodalLikeOut'
            % Constrain prior to be bimodal
            numComps  = randi(maxComps-1) + 1;
    
            priorWs   = rand(1,numComps);
            priorWs   = priorWs/sum(priorWs);
            
            if numComps < 4
                priorMus  = [rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
                             rand(1,numComps-1)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
            else
                priorMus  = [rand(1,2)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
                             rand(1,numComps-2)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
            end
            priorSigs = max([sig1 sig2])*(1 + rand(1,numComps)*prScale*0.7);
            
            % Reselect stimulus means, such that they don't fall between modes of prior
            choseOne = randi(2,[1 2]);
            
            leftBound  = min(priorMus);
            rightBound = max(priorMus);
            leftBox    = abs(stimBnds(1)-leftBound);
            rightBox   = abs(stimBnds(2)-rightBound);
            
            pull1      = [rand*leftBox - stimBnds(1), rand*rightBox - rightBound];
            pull2      = [rand*leftBox - stimBnds(1), rand*rightBox - rightBound];
            
            x1   = pull1(choseOne(1));
            x2   = pull2(choseOne(2));
            
            mu1  = x1;
            mu2  = x2;
            
        case 'none'
            % Unconstrained prior
            numComps  = randi(maxComps-1) + 1;
    
            priorWs   = rand(1,numComps);
            priorWs   = priorWs/sum(priorWs);
            priorMus  = rand(1,numComps)*stimBnds(2) - stimBnds(2)/2;
            priorSigs = rand(1,numComps)*prScale;
            
    end

    % Define support centered on 2D likelihood peak, extend to 4SD (99.994% AUC) and
    % choose granularity of 8SD units/dx
    dx       = 300;
    
    suppLB1  = x1-4*sig1;
    suppUB1  = x1+4*sig1;
    support1 = linspace(suppLB1,suppUB1,dx);
    
    suppLB2  = x2-4*sig2;
    suppUB2  = x2+4*sig2;
    support2 = linspace(suppLB2,suppUB2,dx);
    
    % Make the prior and determine the best-fitting single Gaussian
    suppLBP     = min(priorMus) - 3*max(priorSigs);
    suppUBP     = max(priorMus) + 3*max(priorSigs);
    suppPrior   = linspace(suppLBP,suppUBP,dx);
    
    [thisPrior]  = getMoGPrior(priorSigs,priorMus,priorWs,suppPrior);
    [sgMu,sgSig] = fitSingleGauss(suppPrior,thisPrior);
    
    %%%%%%% Debug
    if mod(ii,numSims/5) == 0
        figure;
        hold on;
        plot(suppPrior,thisPrior);
        sgfit = normpdf(suppPrior,sgMu,sgSig);
        plot(suppPrior,sgfit/sum(sgfit));
        set(gca,'yscale','log');
    end
    %%%%%%%%%%%%%%%%
    
    
    % Numerically calculate psychometric function
    [pFxnNum] = calcMoGPFxn_Numeric(support1,support2,priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2,0);
    
    % Analytically calculate psychometric function
    [pFxnAna] = calcMoGPFxn_Analytic(priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2);
    
    % Analytically calculate psychometric function from best-fitting single Gaussian
    % (analytical approximation reduces down to exact solution to single
    % Gaussian for a single component)
    [pFxnSG]  = calcMoGPFxn_Analytic(sgMu,sgSig,[1],mu1,sig1,mu2,sig2);
    
    % Calculate raw error
    error   = pFxnAna - pFxnNum;
    errorSG = pFxnSG - pFxnNum;
    
    % Package this round of the simulation up
    simMat(ii).x1        = x1;
    simMat(ii).mu1       = mu1;
    simMat(ii).sig1      = sig1;
    simMat(ii).x2        = x2;
    simMat(ii).mu2       = mu2;
    simMat(ii).sig2      = sig2;
    simMat(ii).priorWs   = priorWs;
    simMat(ii).priorMus  = priorMus;
    simMat(ii).priorSigs = priorSigs;
    simMat(ii).sgMu      = sgMu;
    simMat(ii).sgSig     = sgSig;
    simMat(ii).pFxnAna   = pFxnAna;
    simMat(ii).pFxnNum   = pFxnNum;
    simMat(ii).pFxnSG    = pFxnSG;
    simMat(ii).error     = error;
    simMat(ii).errorSG   = errorSG;
    
end

%% Histogram results

% Grab errors and distances between stimuli
stimDists = arrayfun(@(x) x.x2 - x.x1,simMat);
estErrors = arrayfun(@(x) x.error,simMat);
estErrSG  = arrayfun(@(x) x.errorSG,simMat);

% Histogram errors
% maxError  = max(abs([estErrors estErrSG]));
% maxError  = 0.3;
maxError  = 0.05;
binEdges  = linspace(-maxError,maxError,50);
rmsError  = sqrt(nanmean(estErrors.^2));
rmsErrSG  = sqrt(nanmean(estErrSG.^2));


%% Plot results

% Plot two approximations 
f1          = figure;
f1.Position = [100 100 650 650];
hold on;

pAna = arrayfun(@(x) x.pFxnAna,simMat);
pNum = arrayfun(@(x) x.pFxnNum,simMat);
pSG  = arrayfun(@(x) x.pFxnSG,simMat);

plot([0 1],[0 1],'--k','linewidth',4);
s1 = scatter(pNum,pSG,60,colorMat(2,:),'filled');
s2 = scatter(pNum,pAna,60,colorMat(1,:),'filled');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1],'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[0 0.5 1]);
xlabel('Numerical approximation');
ylabel('Analytical approximation');
legend([s2,s1],{'Analytical','Single gaussian'},'location','northwest');

% Error histogram
f2          = figure;
f2.Position = [500 100 750 700];
hold on;

hist1       = histogram(estErrors,'BinEdges',binEdges,'facecolor',colorMat(1,:),'edgecolor','none');
hist2       = histogram(estErrSG,'BinEdges',binEdges,'facecolor',colorMat(2,:),'edgecolor','none');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[binEdges(1) binEdges(end)],'ylim',...
            [0 max([hist1.Values hist2.Values])*1.05]);
xlabel('Analytical - Numerical approximation ');
ylabel('Count');
text(binEdges(1)*0.88*[1 1],max([hist1.Values hist2.Values])*0.88*[1 1],['RMS error = ',...
    num2str(round(rmsError,2,'significant'))],'fontsize',15);
text(binEdges(1)*0.88*[1 1],max([hist1.Values hist2.Values])*0.83*[1 1],['RMS error (SG) = ',...
    num2str(round(rmsErrSG,2,'significant'))],'fontsize',15);
legend([hist1,hist2],{'Analytical','Single gaussian'},'location','northeast');
