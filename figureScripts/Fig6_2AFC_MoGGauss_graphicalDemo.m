%% Figure 6: 2AFC Demo - MoG Prior

% Note: this script take a few minutes to run due to numerically
% calculating the psychometric function values

clear all
close all

% Define numerical support
suppLB  = 0;                                % lower bound
suppUB  = 5;                                % upper bound
numBins = 300;                              % number of bins
supp    = linspace(suppLB,suppUB,numBins);  % stimulus support
dx      = diff(supp(1:2));                  % precision

% Define observer/stimulus parameters
%------------------% 
% Stimulus values
x1 = 3;
x2 = 3;

% Prior
pGam = [2;0.6];
pNu  = [0;0];
pW   = [0.5;0.5];

prior  = getMoGPrior(pGam,pNu,pW,supp);
%------------------% 

% Toggle on/off saving figures
saveOn = 0;


%% Stimulus Combination 1

% Define stimulus standard deviations
s1a = 0.75;
s2a = 0.5;

% Generate likelihoods
p2a = normpdf(supp,x1,s1a);
p3a = normpdf(supp,x2,s2a);

% Numerically calculate probability observer reports stimulus 2 larger
% p(yes|x1,x2). Note, the last argument tells the function to plot a figure
% of the joint likelihood function and decision boundary.
[pf1,~,f1] = get2AFCresponseP(x1,x2,[s1a s2a],pW,pNu,pGam,1);


%% Stimulus Combination 2
s1b = 0.637;
s2b = 0.637;

p2b = normpdf(supp,x1,s1b);
p3b = normpdf(supp,x2,s2b);

[pf2,~,f2] = get2AFCresponseP(x1,x2,[s1b s2b],pW,pNu,pGam,1);


%% Stimulus Combination 3
s1c = 0.5;
s2c = 0.75;

p2c = normpdf(supp,x1,s1c);
p3c = normpdf(supp,x2,s2c);

[pf3,~,f3] = get2AFCresponseP(x1,x2,[s1c s2c],pW,pNu,pGam,1);


%% Psychometric Functions

% Define numerical support (values of stimulus 2) for psychometric curves
numX2 = 15;
twoAFCStim = linspace(0,6,numX2);

% Initialize psychometric curve 
pf1_plot = nan(numX2,1);
pf2_plot = nan(numX2,1);
pf3_plot = nan(numX2,1);

% Loop over different values of stimulus 2 to make a "full" curve
for ii = 1:numel(twoAFCStim)
    
    if ii == (numX2/2 + 0.5)
        pf1_plot(ii) = pf1;
        pf2_plot(ii) = pf2;
        pf3_plot(ii) = pf3;
    else
        pf1_plot(ii) = get2AFCresponseP(x1,twoAFCStim(ii),[s1a s2a],pW,pNu,pGam,0);
        pf2_plot(ii) = get2AFCresponseP(x1,twoAFCStim(ii),[s1b s2b],pW,pNu,pGam,0);
        pf3_plot(ii) = get2AFCresponseP(x1,twoAFCStim(ii),[s1c s2c],pW,pNu,pGam,0);
    end
    
end


%% Plot figures

% Likelihood combination 1
%------------------% 
f1.Position = [100 100 650 650];

f1b = figure;
f1b.Position = [100 200 650 325];
hold on;

pl1a = plot(supp,prior,'r','linewidth',4);
pl2a = plot(supp,p2a,'k','linewidth',4);
pl3a = plot(supp,p3a,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 max([prior(:); p2a(:); p3a(:)])]);
legend([pl1a pl2a pl3a],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');


% Likelihood combination 2
%------------------% 
f2.Position = [300 100 650 650];

f2b = figure;
f2b.Position = [300 200 650 325];
hold on;

pl1b = plot(supp,prior,'r','linewidth',4);
pl2b = plot(supp,p2b,'k','linewidth',4);
pl3b = plot(supp,p3b-0.01,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 max([prior(:); p2b(:); p3b(:)])]);
legend([pl1b pl2b pl3b],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');


% Likelihood combination 3
%------------------% 
f3.Position = [500 100 650 650];

f3b = figure;
f3b.Position = [500 200 650 325];
hold on;

pl1c = plot(supp,prior,'r','linewidth',4);
pl2c = plot(supp,p2c,'k','linewidth',4);
pl3c = plot(supp,p3c,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',20,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 max([prior(:); p2c(:); p3c(:)])]);
legend([pl1c pl2c pl3c],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');


% Psychometric function
%------------------% 
f4 = figure;
f4.Position = [400 100 650 650];
hold on;

plot(twoAFCStim,pf1_plot,'k','linewidth',4);
plot(twoAFCStim,pf2_plot,'k','linewidth',4);
plot(twoAFCStim,pf3_plot,'k','linewidth',4);
scatter(x2,pf1,200,'k','filled');
scatter(x2,pf2,200,'k','filled');
scatter(x2,pf3,200,'k','filled');

set(gca,'fontsize',40,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[suppLB 6],...
        'ytick',[0 0.5 1],'ylim',[0 1]);
xlabel('Stimulus 2 (x_{2})');
ylabel('p(yes|x_{1},x_{2})');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig6_2AFC_MoGGauss_graphicalDemo'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig6/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'like2Da.svg']);
    saveas(f1b,[sDir,'like1Da.svg']);
    saveas(f2,[sDir,'like2Db.svg']);
    saveas(f2b,[sDir,'like1Db.svg']);
    saveas(f3,[sDir,'like2Dc.svg']);
    saveas(f3b,[sDir,'like1Dc.svg']);
    saveas(f4,[sDir,'psychFxns.svg']);
    
end