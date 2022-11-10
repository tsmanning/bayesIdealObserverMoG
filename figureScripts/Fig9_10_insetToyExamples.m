%% Fig 9/10 - inset of prior/likelihood selection process

% Note: this script randomizes the stimuli and prior so the plots will
% differ from those in the paper and each previous run. 

clear all;
close all;

% Define color palette used for plotting
colorMat = colororder;

%% Setup simulation pars

% Define bounds for stimuli/measurements
stimBnds    = [-1 1];

% Set prior scale factor (e.g. 1 - prior components are 1-2x wider than
% likelihoods, 2 - prior components are 1-4x wider, etc)
prScale     = 3;

% Set max number of components in prior
% maxComps    = 5;
maxComps    = 2;


%% Run simulation

% Randomly select stimuli, define random measurement distributions
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

suppLB2  = x2-4*sig2;
suppUB2  = x2+4*sig2;
support2 = linspace(suppLB2,suppUB2,dx);

% Define a random prior
% Constrain prior s.t. it has a long tail and means are both zero
numComps1  = randi(maxComps-1) + 1;

priorWs1   = rand(1,numComps1);
priorWs1   = priorWs1/sum(priorWs1);
priorMus1  = zeros(1,numComps1);
priorSigs1 = [1.1*max([sig1 sig2]) max([sig1 sig2])*(1.75 + rand(1,numComps1-1)*prScale)];

% Constrain prior to be bimodal
numComps2  = randi(maxComps-1) + 1;

priorWs2   = 0.3*ones(1,numComps2)+rand(1,numComps2);
priorWs2   = priorWs2/sum(priorWs2);

priorSigs2 = max([sig1 sig2])*(1 + rand(1,numComps2)*prScale*0.125);

if numComps2 < 4
    %                 priorMus  = [rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
    %                              rand(1,numComps-1)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
    
    priorMus2(1)  = rand(1,1)*stimBnds(1)*0.5 + stimBnds(1)*0.5;
    priorMus2(2)  = priorMus2(1) + 3*max(priorSigs2) + rand(1)*stimBnds(2)*0.4;
else
    priorMus2  = [rand(1,2)*stimBnds(1)*0.5 + stimBnds(1)*0.5,...
        rand(1,numComps2-2)*stimBnds(2)*0.5 + stimBnds(2)*0.5];
end

% Make the prior and determine the best-fitting single Gaussian
suppLBP1     = min(priorMus1) - 3*max(priorSigs1);
suppUBP1     = max(priorMus1) + 3*max(priorSigs1);
suppPrior1   = linspace(suppLBP1,suppUBP1,dx);

[thisPrior1]  = getMoGPrior(priorSigs1,priorMus1,priorWs1,suppPrior1);

suppLBP2     = min(priorMus2) - 3*max(priorSigs2);
suppUBP2     = max(priorMus2) + 3*max(priorSigs2);
suppPrior2   = linspace(suppLBP2,suppUBP2,dx);

[thisPrior2]  = getMoGPrior(priorSigs2,priorMus2,priorWs2,suppPrior2);


%% Plot
% Plot priors and likelihoods
f1 = figure;
f1.Position = [100 100 650 650];
hold on

like1 = normpdf(suppPrior1,mu1,sig1);
like2 = normpdf(suppPrior1,mu2,sig2);

p1 = plot(suppPrior1,thisPrior1/(sum(thisPrior1)*dx),'r','linewidth',4);
p2 = plot(suppPrior1,like1/(sum(like1)*dx),'k','linewidth',4);
p3 = plot(suppPrior1,like2/(sum(like2)*dx),'--k','linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[suppLBP1,suppUBP1],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([p1 p2 p3],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest');

% Plot priors and likelihoods
f2 = figure;
f2.Position = [700 100 650 650];
hold on

like1 = normpdf(suppPrior2,mu1,sig1);
like2 = normpdf(suppPrior2,mu2,sig2);

p1 = plot(suppPrior2,thisPrior2/(sum(thisPrior2)*dx),'r','linewidth',4);
p2 = plot(suppPrior2,like1/(sum(like1)*dx),'k','linewidth',4);
p3 = plot(suppPrior2,like2/(sum(like2)*dx),'--k','linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[suppLBP2,suppUBP2],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([p1 p2 p3],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest');

