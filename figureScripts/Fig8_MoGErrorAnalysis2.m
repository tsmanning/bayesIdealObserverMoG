%% Figure 8 - Error in 2AFC approximation: bimodal prior

clear all
close all

% Define color palette used for plotting
colorMat = colororder;

% Toggle on/off saving figures
saveOn = 0;


%% Load in data from error analysis

% Run testAnalyError.m (with priorConsts set to 'bimodal') to obtain this data
simDir = '/media/tyler/Data/MATLAB/cooperLab/4-Papers/MoGpaper/';
load([simDir,'bimodalPrior_noConstraints_new'],'estErrors','estErrSG','simMat')


%% Error scatters

% Extract estimates of psychometric function from each run of the
% simulation for the MoG analytical approximation, numerical calculation, and
% single gaussian approximation
pAna = arrayfun(@(x) x.pFxnAna,simMat);
pNum = arrayfun(@(x) x.pFxnNum,simMat);
pSG  = arrayfun(@(x) x.pFxnSG,simMat);


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

% Define bounds for stimuli/measurements
stimBnds    = [-1 1];

% Set prior scale factor (e.g. 1 - prior components are 1-2x wider than
% likelihoods, 2 - prior components are 1-4x wider, etc)
prScale     = 3;

% Set max number of components in prior
maxComps    = 2;

% Randomly select two stimuli, define random measurement distributions
x1   = rand*stimBnds(2)*2 - stimBnds(2);
mu1  = x1;
sig1 = rand*stimBnds(2)*0.5;

x2   = rand*stimBnds(2)*2 - stimBnds(2);
mu2  = x2;
sig2 = rand*stimBnds(2)*0.5;

% Define support centered on 2D likelihood peak, extend to 4SD (99.994% AUC) and
% choose granularity of 8SD units/dx
dx       = 300;

suppLB2  = x2-4*sig2;
suppUB2  = x2+4*sig2;
support2 = linspace(suppLB2,suppUB2,dx);

% Define a random prior

% Constrain prior to be bimodal
numComps2  = randi(maxComps-1) + 1;

priorWs2   = 0.3*ones(1,numComps2)+rand(1,numComps2);
priorWs2   = priorWs2/sum(priorWs2);

priorSigs2 = max([sig1 sig2])*(1 + rand(1,numComps2)*prScale*0.125);

if numComps2 < 4
    priorMus2(1)  = rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5;
    priorMus2(2)  = priorMus2(1) + 3*max(priorSigs2) + rand(1)*stimBnds(2)*0.4;
else
    priorMus2  = [rand(1,2)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
        rand(1,numComps2-2)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
end

% Make the prior
suppLBP2     = min(priorMus2) - 3*max(priorSigs2);
suppUBP2     = max(priorMus2) + 3*max(priorSigs2);
suppPrior2   = linspace(suppLBP2,suppUBP2,dx);

[thisPrior2]  = getMoGPrior(priorSigs2,priorMus2,priorWs2,suppPrior2);


%% Plot

% Error scatter single Gaussian fit
%-----------------------------%
f1 = figure;
f1.Position = [100 300 600 500];
hold on;

plot([0 1],[0 1],'--k','linewidth',4);
s1 = scatter(pNum,pSG,60,colorMat(2,:),'filled');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1],'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[0 0.5 1]);
xlabel('Numerical evaluation');
ylabel('SG approximation');
legend([s1],{'Single Gaussian'},'location','northwest','box','off');

% Error scatter mixture of Gaussian fit
%-----------------------------%
f2 = figure;
f2.Position = [100 300 600 500];
hold on;

plot([0 1],[0 1],'--k','linewidth',4);
s2 = scatter(pNum,pAna,60,colorMat(1,:),'filled');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1],'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[0 0.5 1]);
xlabel('Numerical evaluation');
ylabel('MoG approximation');
legend([s2],{'MoG approximation'},'location','northwest','box','off');

% Mean signed error 
%-----------------------------%
f3 = figure;
f3.Position = [100 100 850 750];
hold on;

plot([0 1],[0 0],'color',[0.15 0.15 0.15]);
p1 = plot(binCents,MSE,'color',colorMat(1,:),'linewidth',5);
p2 = plot(binCents,MSE_SG,'color',colorMat(2,:),'linewidth',5);
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
p4 = plot(binCents,RMSE_SG,'color',colorMat(2,:),'linewidth',5);
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
    
    splPath = regexp(which('Fig8_MoGErrorAnalysis2'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig8/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'bimodalErrorScatterSG.svg']);
    saveas(f2,[sDir,'bimodalErrorScatterMOG.svg']);
    saveas(f3,[sDir,'bimodalMSE.svg']);
    saveas(f4,[sDir,'bimodalRMSE.svg']);
    saveas(f5,[sDir,'bimodalToyExample.svg']);
    
end
