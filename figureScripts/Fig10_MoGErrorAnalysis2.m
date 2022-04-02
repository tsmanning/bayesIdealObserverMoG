%% Figure 10 - Error in 2AFC approximation: bimodal prior

clear all
close all

% Define color palette used for plotting
colorMat = colororder;

% Toggle on/off saving figures
saveOn = 0;


%% Load in data from error analysis

% Run testAnalyError.m (with priorConsts set to 'bimodal') to obtain this data
splPath = regexp(which('Fig9_MoGErrorAnalysis'),filesep,'split');
topDir  = [fullfile(splPath{1:numel(splPath)-2}),filesep];

load([topDir,'bimodalPrior'],'simMat')


%% Error scatters

% Extract estimates of psychometric function from each run of the
% simulation for the MoG analytical approximation, numerical calculation, and
% single gaussian approximation
pAna = arrayfun(@(x) x.pFxnAna,simMat);
pNum = arrayfun(@(x) x.pFxnNum,simMat);
pSG  = arrayfun(@(x) x.pFxnSG,simMat);

estErrors = arrayfun(@(x) x.error,simMat);
estErrSG  = arrayfun(@(x) x.errorSG,simMat);


%% Calculate signed and unsigned error for the analytical approximation

% Define bin edges over which we will histogram errors (values of stimulus 2)
numBins   = 20;
pFxnEdges = linspace(0,1,numBins + 1);
binCents  = pFxnEdges(1:end-1) + diff(pFxnEdges(1:2))/2;

% Get the mean signed error and RMS error for each bin
MSE = nan(numBins,1);
MSE_SG = nan(numBins,1);
RMSE = nan(numBins,1);
RMSE_SG = nan(numBins,1);

for ii = 1:numBins

    theseInds  = (pNum>pFxnEdges(ii)) & (pNum<pFxnEdges(ii+1));
    MSE(ii)    = mean(estErrors(theseInds));
    MSE_SG(ii) = nanmean(estErrSG(theseInds));
    
    RMSE(ii)    = sqrt(mean(estErrors(theseInds).^2));
    RMSE_SG(ii) = sqrt(nanmean(estErrSG(theseInds).^2));
    
end

peakMSE  = round(1.4*max(abs([MSE MSE_SG])),1,'significant');
% peakRMSE = round(1.4*max(abs([RMSE RMSE_SG])),1,'significant');
peakRMSE = 0.12;


%% Generate an example prior and likelihood pair like those used in error simulation

% Select an index from the error simulation
exind = 83;

% Define bounds for stimuli/measurements
stimBnds    = [-1 1];

% Set prior scale factor (e.g. 1 - prior components are 1-2x wider than
% likelihoods, 2 - prior components are 1-4x wider, etc)
prScale     = 3;

% Set max number of components in prior
maxComps    = 2;

% Get noise and stimulus parameters for this simulation run
mu1  = simMat(exind).mu1;
sig1 = simMat(exind).sig1;

mu2  = simMat(exind).mu2;
sig2 = simMat(exind).sig2;

% Define support granularity
dx       = 300;

% Get prior parameters for this simulation run
priorWs2   = simMat(exind).priorWs;
priorMus2  = simMat(exind).priorMus;
priorSigs2 = simMat(exind).priorSigs;

% Make the prior
suppLBP2     = min(priorMus2) - 3*max(priorSigs2);
suppUBP2     = max(priorMus2) + 3*max(priorSigs2);
suppPrior2   = linspace(suppLBP2,suppUBP2,dx);

[thisPrior2]  = getMoGPrior(priorSigs2,priorMus2,priorWs2,suppPrior2);


%% Fill out res of psychometric function for other test velocities

numTestVels = 20;

pFxnNum = nan(numTestVels,1);
pFxnAna = nan(numTestVels,1);
pFxnSG  = nan(numTestVels,1);
    
exSig1    = simMat(exind).sig1;
exSig2    = simMat(exind).sig2;

priorMus  = simMat(exind).priorMus;
priorSigs = simMat(exind).priorSigs;
priorWs   = simMat(exind).priorWs;

mu1       = simMat(exind).mu1;

[~,maxComp] = max(simMat(exind).priorWs);
maxSig      = max([exSig1 exSig2]);
pseEst      = mu1*(simMat(exind).priorSigs(maxComp)^2 + exSig2^2)/(simMat(exind).priorSigs(maxComp)^2 + exSig1^2);
testMus     = linspace(pseEst - 4*max([exSig1 exSig2]), pseEst + 4*max([exSig1 exSig2]),numTestVels-1);
testMus     = sort([testMus simMat(exind).x2]);

suppLB1  = mu1 - 4*exSig1;
suppUB1  = mu1 + 4*exSig1;
support1 = linspace(suppLB1,suppUB1,dx);

sgMu     = simMat(exind).sgMu;
sgSig    = simMat(exind).sgSig;

for ii = 1:numTestVels
    
    thisTestMu = testMus(ii);
    
    suppLB2  = thisTestMu - 4*exSig2;
    suppUB2  = thisTestMu + 4*exSig2;
    support2 = linspace(suppLB2,suppUB2,dx);
    
    % Numerically calculate psychometric function
    pFxnNum(ii) = calcMoGPFxn_Numeric(support1,support2,priorMus,priorSigs,priorWs,mu1,exSig1,thisTestMu,exSig2,0);
    
    % Analytically calculate psychometric function
    pFxnAna(ii) = calcMoGPFxn_Analytic(priorMus,priorSigs,priorWs,mu1,exSig1,thisTestMu,exSig2);
    
    % Analytically calculate psychometric function from best-fitting single Gaussian
    % (analytical approximation reduces down to exact solution to single
    % Gaussian for a single component)
    pFxnSG(ii)  = calcMoGPFxn_Analytic(sgMu,sgSig,[1],mu1,exSig1,thisTestMu,exSig2);
    
end


%% Plot

% Scatter with numerical evaluation/MoG/SG on same plot
%-----------------------------%
f0 = figure;
f0.Position = [200 800 600 500];
hold on;

plot(testMus,pFxnNum,'k','linewidth',4);
plot(testMus,pFxnAna,'color',colorMat(1,:),'linewidth',4);
plot(testMus,pFxnSG,'color',colorMat(3,:),'linewidth',4);

scatter(simMat(exind).x2,simMat(exind).pFxnNum,120,[0 0 0],'filled');
scatter(simMat(exind).x2,simMat(exind).pFxnAna,120,colorMat(1,:),'filled');
scatter(simMat(exind).x2,simMat(exind).pFxnSG,120,colorMat(3,:),'filled');

plot([testMus(1) simMat(exind).x2],simMat(exind).pFxnNum*[1 1],'--k','linewidth',2);
plot([testMus(1) simMat(exind).x2],simMat(exind).pFxnAna*[1 1],'color',colorMat(1,:),'linewidth',2,'linestyle','--');
plot([testMus(1) simMat(exind).x2],simMat(exind).pFxnSG*[1 1],'color',colorMat(3,:),'linewidth',2,'linestyle','--');

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[testMus(1) testMus(end)],'ylim',[0 1],'ytick',[0 0.5 1]);
xlabel('Stimulus (x_{2})');
ylabel('p("x_{2} > x_{1}")');

% Error scatter single Gaussian fit
%-----------------------------%
f1 = figure;
f1.Position = [100 300 600 500];
hold on;

plot([0 1],[0 1],'--k','linewidth',4);
s1 = scatter(pNum,pSG,60,colorMat(3,:),'filled');
scatter(pNum(exind),pSG(exind),120,[0 0 0],'linewidth',4);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1],'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[0 0.5 1]);
xlabel('Numerical evaluation');
ylabel('SG approximation');
% legend([s1],{'Single Gaussian'},'location','northwest','box','off');

% Error scatter mixture of Gaussian fit
%-----------------------------%
f2 = figure;
f2.Position = [100 300 600 500];
hold on;

plot([0 1],[0 1],'--k','linewidth',4);
s2 = scatter(pNum,pAna,60,colorMat(1,:),'filled');
scatter(pNum(exind),pAna(exind),120,[0 0 0],'linewidth',4);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1],'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[0 0.5 1]);
xlabel('Numerical evaluation');
ylabel('MoG approximation');
% legend([s2],{'MoG approximation'},'location','northwest','box','off');

% Mean signed error 
%-----------------------------%
f3 = figure;
f3.Position = [100 100 850 750];
hold on;

plot([0 1],[0 0],'color',[0.15 0.15 0.15]);
p1 = plot(binCents,MSE,'color',colorMat(1,:),'linewidth',5);
p2 = plot(binCents,MSE_SG,'color',colorMat(3,:),'linewidth',5);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'ylim',0.12*[-1 1],'xtick',[0 0.5 1],'ytick',0.12*linspace(-1,1,5));
xlabel('Numerical evaluation');
ylabel('Mean signed error');
legend([p1 p2],{'MoG approximation','Single Gaussian'},'location','northwest','box','off');

% Unsigned (RMS) error 
%-----------------------------%
f4 = figure;
f4.Position = [1000 100 850 750];
hold on;

plot([0 1],[0 0],'color',[0.15 0.15 0.15]);
p3 = plot(binCents,RMSE,'color',colorMat(1,:),'linewidth',5);
p4 = plot(binCents,RMSE_SG,'color',colorMat(3,:),'linewidth',5);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'ylim',[0 0.12],'xtick',[0 0.5 1],'ytick',0.12*linspace(0,1,5));
xlabel('Numerical evaluation');
ylabel('RMS error');
legend([p3 p4],{'MoG approximation','Single Gaussian'},'location','northwest','box','off');

% 2AFC Demo with example prior
%-----------------------------%
f5 = figure;
f5.Position = [100 100 650 650];
hold on;

like1 = normpdf(suppPrior2,mu1,sig1);
like2 = normpdf(suppPrior2,mu2,sig2);

p1 = plot(suppPrior2,thisPrior2/(sum(thisPrior2)*dx),'r','linewidth',4);
p2 = plot(suppPrior2,like1/(sum(like1)*dx),'k','linewidth',4);
p3 = plot(suppPrior2,like2/(sum(like2)*dx),'--k','linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[suppLBP2,suppUBP2],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([p1 p2 p3],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest');


%% Save figures

if saveOn

    sDir = [topDir,'figuresImgs/fig10/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'bimodalErrorScatterSG.svg']);
    saveas(f2,[sDir,'bimodalErrorScatterMOG.svg']);
    saveas(f3,[sDir,'bimodalMSE.svg']);
    saveas(f4,[sDir,'bimodalRMSE.svg']);
    saveas(f5,[sDir,'bimodalToyExample.svg']);
    
end
