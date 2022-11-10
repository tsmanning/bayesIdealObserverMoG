%% Figure 4: 2AFC demo - Single Gaussian Prior

clear all
close all

% Define numerical support
suppLB  = 0;                                % lower bound
suppUB  = 5;                                % upper bound
numBins = 300;                              % number of bins
supp    = linspace(suppLB,suppUB,numBins);  % stimulus support
dx      = diff(supp(1:2));                  % precision

% Define stimulus values
x1 = 3;
x2 = 3;

% Define prior
pGam = 1.5;     % Standard deviation
pNu  = 0;       % Mean

% Save figure images?
saveOn = 0;


%% Stimulus Combination 1

% Define stimulus standard deviations
s1a = 0.75;     
s2a = 0.5;

% Define shrinkage factors for the two stimuli
alpha1a = (pGam^2)/(s1a^2 + pGam^2);
alpha2a = (pGam^2)/(s2a^2 + pGam^2);

% Create a joint likelihood distribution to plot
like2Da = normpdf(supp,x2,s2a)'*normpdf(supp,x1,s1a);

% Create the single Gaussian prior and the two likelihoods
prior = normpdf(supp,pNu,pGam);
like1a = normpdf(supp,x1,s1a);
like2a = normpdf(supp,x2,s2a);


%% Combination 2
s1b = 0.637;
s2b = 0.637;

alpha1b = (pGam^2)/(s1b^2 + pGam^2);
alpha2b = (pGam^2)/(s2b^2 + pGam^2);

like2Db = normpdf(supp,x2,s2b)'*normpdf(supp,x1,s1b);

like1b = normpdf(supp,x1,s1b);
like2b = normpdf(supp,x2,s2b);


%% Combination 3
s1c = 0.5;
s2c = 0.75;

alpha1c = (pGam^2)/(s1c^2 + pGam^2);
alpha2c = (pGam^2)/(s2c^2 + pGam^2);

like2Dc = normpdf(supp,x2,s2c)'*normpdf(supp,x1,s1c);

like1c = normpdf(supp,x1,s1c);
like2c = normpdf(supp,x2,s2c);


%% Psychometric Functions

% Define psychometric function as the cumulative of the difference of two Gaussian random
% variables
psychFxn = @(a1,a2,s1,s2,sig1,sig2) normcdf( (a2*s2 - a1*s1)/sqrt((a2^2)*(sig2^2) + (a1^2)*(sig1^2)) );

% Define support for psychometric function (values of stimulus 2)
twoAFCStim = suppLB:0.1:6;


%% Plot figures

% Define gamma'd colormap to use for plots
cmap = repmat(linspace(0,1,256)'.^1.0,[1 3]);

% Likelihood combination 1
%------------------% 
yMax = max(like2a(:));

f1 = figure;
f1.Position = [100 100 650 650];
hold on;

imagesc(supp,supp,like2Da); axis image; axis xy;
plot([suppLB suppUB],alpha1a/alpha2a*[suppLB suppUB],'w','linewidth',5);
text([0.25],[1.5],'Yes','color',[1 1 1],'fontsize',40);
text([1.25],[0.5],'No','color',[1 1 1],'fontsize',40);
set(gca,'fontsize',40,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[suppLB suppUB]);
xlabel('Stimulus 1 (x_{1})');
ylabel('Stimulus 2 (x_{2})');
colormap(cmap);

f1b = figure;
f1b.Position = [100 200 650 325];
hold on;

pl1a = plot(supp,prior,'r','linewidth',4);
pl2a = plot(supp,like1a,'k','linewidth',4);
pl3a = plot(supp,like2a,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',40,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 yMax]);
legend([pl1a pl2a pl3a],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');

% Likelihood combination 2
%------------------% 
f2 = figure;
f2.Position = [300 100 650 650];
hold on;

imagesc(supp,supp,like2Db); axis image; axis xy;
plot([suppLB suppUB],alpha1b/alpha2b*[suppLB suppUB],'w','linewidth',5);
text([0.25],[1.5],'Yes','color',[1 1 1],'fontsize',40);
text([1.25],[0.5],'No','color',[1 1 1],'fontsize',40);
set(gca,'fontsize',40,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[suppLB suppUB]);
xlabel('Stimulus 1 (x_{1})');
ylabel('Stimulus 2 (x_{2})');
colormap(cmap);

f2b = figure;
f2b.Position = [300 200 650 325];
hold on;

pl1b = plot(supp,prior,'r','linewidth',4);
pl2b = plot(supp,like1b+0.015,'k','linewidth',4);
pl3b = plot(supp,like2b,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',40,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 yMax]);
legend([pl1b pl2b pl3b],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');

% Likelihood combination 3
%------------------% 
f3 = figure;
f3.Position = [300 100 650 650];
hold on;

imagesc(supp,supp,like2Dc); axis image; axis xy;
plot([suppLB suppUB],alpha1c/alpha2c*[suppLB suppUB],'w','linewidth',5);
text([0.25],[1.5],'Yes','color',[1 1 1],'fontsize',40);
text([1.25],[0.5],'No','color',[1 1 1],'fontsize',40);
set(gca,'fontsize',40,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[suppLB suppUB]);
xlabel('Stimulus 1 (x_{1})');
ylabel('Stimulus 2 (x_{2})');
colormap(cmap);

f3b = figure;
f3b.Position = [300 200 650 325];
hold on;

pl1c = plot(supp,prior,'r','linewidth',4);
pl2c = plot(supp,like1c,'k','linewidth',4);
pl3c = plot(supp,like2c,'color',[0.8 0.8 0.8],'linewidth',4);

set(gca,'fontsize',40,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[suppLB suppUB],...
        'ytick',[],'ylim',[0 yMax]);
legend([pl1c pl2c pl3c],{'Prior','Likelihood 1','Likelihood 2'},'location','northwest','box','off');

% Psychometric function
%------------------% 
f4 = figure;
f4.Position = [400 100 650 650];
hold on;

plot(twoAFCStim,psychFxn(alpha1a,alpha2a,x1,twoAFCStim,s1a,s2a),'k','linewidth',4);
plot(twoAFCStim,psychFxn(alpha1b,alpha2b,x1,twoAFCStim,s1b,s2b),'k','linewidth',4);
plot(twoAFCStim,psychFxn(alpha1c,alpha2c,x1,twoAFCStim,s1c,s2c),'k','linewidth',4);
scatter(x2,psychFxn(alpha1a,alpha2a,x1,x2,s1a,s2a),200,'k','filled');
scatter(x2,psychFxn(alpha1b,alpha2b,x1,x2,s1b,s2b),200,'k','filled');
scatter(x2,psychFxn(alpha1c,alpha2c,x1,x2,s1c,s2c),200,'k','filled');

set(gca,'fontsize',40,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[suppLB 6],...
        'ytick',[0 0.5 1],'ylim',[0 1]);
xlabel('Stimulus 2 (x_{2})');
ylabel('p(yes|x_{1},x_{2})');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig4_2AFC_SinGauss_graphicalDemo'),filesep,'split');
    topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig4/'];
    
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