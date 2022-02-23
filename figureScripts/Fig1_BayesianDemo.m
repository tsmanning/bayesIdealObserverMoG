%% Figure 1: Non-Gaussian Bayesian computation demo figure

clear all 
close all

% Define numerical support
dx  = 0.001;        % precision of support
x   = -1.5:dx:1.5;  % range of support

% Observer/stim params
%------------------% 
pNu  = [-0.3 0.2 -0.3];     % prior component means
pGam = [0.3 0.4 0.8];       % prior component standard deviations
pW   = [0.3 0.3 0.4];       % prior component weights

lMu  = 0.5;                 % likelihood mean
lSig = 0.3;                 % likelihood standard deviation
%------------------% 

% Toggle on/off saving figures
saveOn = 0;


%% Generate distributions
prior = getMoGPrior(pGam,pNu,pW,x);
like  = normpdf(x,lMu,lSig);

post = like.*prior;         % posterior is product of prior and likelihood
post = post/(sum(post)*dx); % normalize posterior


%% Plot
f = figure;
f.Position = [100 100 1000 500];
hold on;

p1 = plot(x,prior,'color',[0.9 0 0],'linewidth',5);
p2 = plot(x,like,'color',[0 0 0],'linewidth',5);
p3 = plot(x,post,'color',[0 0 1],'linewidth',5);

set(gca,'fontsize',20,'plotboxaspectratio',[2 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([p1,p2,p3],{'Prior, p(x)','Likelihood, p(m|x)','Posterior, p(x|m)'},...
    'location','northwest','fontsize',20);
legend('boxoff');


%% Save figures
if saveOn
    
    splPath = regexp(which('Fig1_BayesianDemo'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir    = [topDir,'figuresImgs/fig1/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f,[sDir,'Fig1_bayesDemo.svg']);
    
end