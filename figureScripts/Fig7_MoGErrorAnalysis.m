%% Figure 7 - Error in 2AFC approximation: zero-centered heavy tailed MoG

clear all
close all

% Define color palette used for plotting
colorMat = colororder;

% Toggle on/off saving figures
saveOn = 0;


%% Load in data from error analysis

% Run testAnalyError.m (with priorConsts set to 'longTail') to obtain this data
simDir = '/media/tyler/Data/MATLAB/cooperLab/4-Papers/MoGpaper/';
load([simDir,'widePrior_wider2ndcomp'],'estErrors','estErrSG','simMat')


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
    MSE_SG(ii) = mean(estErrSG(theseInds));
    
    RMSE(ii)    = sqrt(mean(estErrors(theseInds).^2));
    RMSE_SG(ii) = sqrt(mean(estErrSG(theseInds).^2));
    
end

peakMSE  = round(1.3*max(abs([MSE MSE_SG])),1,'significant');
peakRMSE = round(1.3*max(abs([RMSE RMSE_SG])),1,'significant');


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

suppLB1  = x1-4*sig1;
suppUB1  = x1+4*sig1;
support1 = linspace(suppLB1,suppUB1,dx);

% Define a random prior
% Constrain prior s.t. it has a long tail and means are both zero
numComps1  = randi(maxComps-1) + 1;

priorWs1   = rand(1,numComps1);
priorWs1   = priorWs1/sum(priorWs1);
priorMus1  = zeros(1,numComps1);
priorSigs1 = [1.1*max([sig1 sig2]) max([sig1 sig2])*(1.75 + rand(1,numComps1-1)*prScale)];

% Make the prior and determine the best-fitting single Gaussian
suppLBP1     = min(priorMus1) - 3*max(priorSigs1);
suppUBP1     = max(priorMus1) + 3*max(priorSigs1);
suppPrior1   = linspace(suppLBP1,suppUBP1,dx);

[thisPrior1]  = getMoGPrior(priorSigs1,priorMus1,priorWs1,suppPrior1);


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

like1 = normpdf(suppPrior1,mu1,sig1);
like2 = normpdf(suppPrior1,mu2,sig2);

p1 = plot(suppPrior1,thisPrior1/(sum(thisPrior1)*dx),'r','linewidth',4);
p2 = plot(suppPrior1,like1/(sum(like1)*dx),'k','linewidth',4);
p3 = plot(suppPrior1,like2/(sum(like2)*dx),'--k','linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[suppLBP1,suppUBP1],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([p1 p2 p3],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig7_MoGErrorAnalysis'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig7/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'unimodalErrorScatterSG.svg']);
    saveas(f2,[sDir,'unimodalErrorScatterMOG.svg']);
    saveas(f3,[sDir,'unimodalMSE.svg']);
    saveas(f4,[sDir,'unimodalRMSE.svg']);
    saveas(f5,[sDir,'unimodalToyExample.svg']);
    
end
