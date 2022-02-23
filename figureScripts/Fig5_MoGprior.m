%% Figure 5: MoG prior components to posterior components

clear all
close all

% Toggle on/off saving figures
saveOn = 1;

% Define numerical support
dx = 0.01;
x = -2:dx:2;

% Define observer/stimulus parameters
%------------------% 
% Prior
pNu  = [-0.3 0.2 -0.3];
pGam = [0.3 0.4 0.8];
pW   = [0.3 0.3 0.4];

% Likelihood
lMu  = 0.5;
lSig = 0.3;
%------------------% 

colorMat = colororder;

% Generate full prior
prior = getMoGPrior(pGam,pNu,pW,x);
    
% Generate prior components
p1 = pW(1)*normpdf(x,pNu(1),pGam(1));
p2 = pW(2)*normpdf(x,pNu(2),pGam(2));
p3 = pW(3)*normpdf(x,pNu(3),pGam(3));
    
% Generate likelihood
like  = normpdf(x,lMu,lSig);

% Generate full posterior
post = like.*prior;
post = post/(sum(post)*dx);
    
% Generate posterior components
post1 = like.*p1;
post1 = pW(1)*post1/(sum(post1)*dx); 
post2 = like.*p2;
post2 = pW(2)*post2/(sum(post2)*dx); 
post3 = like.*p3;
post3 = pW(3)*post3/(sum(post3)*dx); 


%% Plot figures

ymax = 1.2*max(like);

f1 = figure;
f1.Position = [100 100 650 650];
hold on;

pr1 = plot(x,p1,'color',[0.8 0.6 0.6],'linewidth',4);
plot(x,p2,'color',[0.8 0.6 0.6],'linewidth',4);
plot(x,p3,'color',[0.8 0.6 0.6],'linewidth',4);
pr = plot(x,prior,'color',[0.8 0 0],'linewidth',4);
li = plot(x,like,'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',30,'xtick',[],'ytick',[],'ylim',[0 ymax]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([pr,pr1,li],{'Prior','Gaussian components','Likelihood'},'location','northwest');

f2 = figure;
f2.Position = [700 100 650 650];
hold on;

po1 = plot(x,post1,'color',[0.6 0.6 1],'linewidth',4);
plot(x,post2,'color',[0.6 0.6 1],'linewidth',4);
plot(x,post3,'color',[0.6 0.6 1],'linewidth',4);
po = plot(x,post,'color',[0 0 1],'linewidth',4);
li = plot(x,like,'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',30,'xtick',[],'ytick',[],'ylim',[0 ymax]);
xlabel('Stimulus (x)');
ylabel('Probability');
legend([po,po1,li],{'Posterior','Gaussian components','Likelihood'},'location','northwest');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig5_MoGprior'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig5/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'MoGpriorComps.svg']);
    saveas(f2,[sDir,'MoGpostComps.svg']);
    
end