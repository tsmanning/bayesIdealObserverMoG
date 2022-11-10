%% Non-Gaussian Bayesian computation demo 2D

clear all
close all

% Define numerical support
dx = 0.001;         % precision
x  = -1.5:dx:1.5;   % stimulus support
m  = x;             % measurement support

gSize = numel(m);   % grid size

% Observer/stim params
%%%%%%%%%%%%%%%%%%
% Prior
pNu  = [-0.3 0.2 -0.3];  % prior component means
pGam = [0.3 0.4 0.8];    % prior component standard deviation 
pW   = [0.3 0.3 0.4];    % prior component weights

% Likelihood
lSig = 0.6;              % likelihood standard deviation
%%%%%%%%%%%%%%%%%%

% Toggle on/off saving figures
saveOn = 0;


%% Define priors

prior   = getMoGPrior(pGam,pNu,pW,x);
    
prior2D = repmat(prior,[numel(m) 1]);


%% Define likelihoods

[xx,mm] = meshgrid(x,m);
like2D  = normpdf(mm,xx,lSig);

% Select two likelihoods and two measurements for plotting
Lind1    = 2401;
Lind2    = 1200;

Mind1    = 2251;
Mind2    = 1501;

like1   = like2D(Lind1,:);
like2   = like2D(Lind2,:);

meas1   = like2D(:,Mind1);
meas2   = like2D(:,Mind2);


%% Define posteriors

post2D = like2D.*prior2D;
post2D = post2D./(repmat(sum(post2D,2),[1 gSize])*dx);

post1  = post2D(Lind1,:);
post2  = post2D(Lind2,:); 



%% Plot everything

% Define gamma'd colormap to use for plots
cmap = repmat(linspace(0,1,256)'.^1.4,[1 3]);

% Prior 1D
f1 = figure;
f1.Position = [100 100 660 275];

plot(x,prior,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');

% Prior 2D
f2 = figure;
f2.Position = [200 100 650 650];

imagesc(x,m,prior2D); axis image; axis xy;
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Measurement (m)');
colormap(cmap);

%---------------------------------%
% Like/Measurement 2D
f3 = figure;
f3.Position = [300 100 650 650];
hold on

imagesc(x,m,like2D); axis image; axis xy;
plot([x(1) x(end)],[x(Lind1) x(Lind1)],'color',[0 0 1],'linewidth',3);
plot([x(1) x(end)],[x(Lind2) x(Lind2)],'color',[0 0 1],'linewidth',3);
plot([x(Mind1) x(Mind1)],[x(1) x(end)],'color',[0 0.8 0],'linewidth',3);
plot([x(Mind2) x(Mind2)],[x(1) x(end)],'color',[0 0.8 0],'linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Measurement (m)');
colormap(cmap);

% Likelihoods 1D
f4a = figure;
f4a.Position = [400 100 660 275];

plot(x,like1,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Likelihood');

f4b = figure;
f4b.Position = [500 100 660 275];

plot(x,like2,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Likelihood');

% Measurements 1D
f5a = figure;
f5a.Position = [600 100 275 650];

plot(x,meas1,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[],'view',[90 -90]);
xlabel('Measurement (m)');
ylabel('Probability');

f5b = figure;
f5b.Position = [700 100 275 650];

plot(x,meas2,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[],'view',[90 -90]);
xlabel('Measurement (m)');
ylabel('Probability');

%---------------------------------%
% Posterior 2D
f6 = figure;
f6.Position = [800 100 650 650];
hold on

imagesc(x,m,post2D); axis image; axis xy;
plot([x(1) x(end)],[x(Lind1) x(Lind1)],'color',[0 0 1],'linewidth',3);
plot([x(1) x(end)],[x(Lind2) x(Lind2)],'color',[0 0 1],'linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Measurement (m)');
colormap(cmap);

% Posteriors 1D
f7a = figure;
f7a.Position = [900 100 660 275];

plot(x,post1,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');

f7b = figure;
f7b.Position = [1000 100 660 275];

plot(x,post2,'k','linewidth',3);
set(gca,'fontsize',20,'plotboxaspectratio',[3 1 1],'xtick',[],'xlim',[x(1) x(end)],'ytick',[]);
xlabel('Stimulus (x)');
ylabel('Probability');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig2_2DBayesianDemo'),filesep,'split');
    topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir    = [topDir,'figuresImgs/fig2/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'prior1D.svg']);
    saveas(f2,[sDir,'prior2D.svg']);
    saveas(f3,[sDir,'like2D.svg']);
    saveas(f4a,[sDir,'like1Da.svg']);
    saveas(f4b,[sDir,'like1Db.svg']);
    saveas(f5a,[sDir,'meas1Da.svg']);
    saveas(f5b,[sDir,'meas1Db.svg']);
    saveas(f6,[sDir,'post2D.svg']);
    saveas(f7a,[sDir,'post1Da.svg']);
    saveas(f7b,[sDir,'post1Db.svg']);
    
end