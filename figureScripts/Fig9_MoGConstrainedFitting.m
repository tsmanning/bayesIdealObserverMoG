%% Figure 9 - Demo showing how to constrain number of parameters

clear all
close all

% Define numerical support
suppLB = -5;
suppUB = 5;
dx     = 300;
supp   = linspace(suppLB,suppUB,dx);

% get color mat
colorMat = colororder;


%% Restrict components to zero mean

% Define number of components to use
numCompsZM = 5;

% Define component parameters
wComp   = rand(1,5);
wComp   = wComp/sum(wComp);
muComp  = zeros(1,numCompsZM);
sigComp = rand(1,5)*suppUB*0.5;

% Make prior from components
priorZM = zeros(1,dx);

for ii = 1:numCompsZM
    
    priorZM = priorZM + wComp(ii)*normpdf(supp,muComp(ii),sigComp(ii));
    
end


%% Restrict components to consistent tiling, fit only the weights of the individual components

% Define number of components to use
numCompsTile  = 10;

% Define centers of each of the components
compCents = linspace(suppLB,suppUB,numCompsTile);

% Calculate component widths based on desired amount of overlap/wiggliness
compWidths = (suppUB-suppLB)/numCompsTile*ones(1,numCompsTile);

% Component weights
% compWeights = rand(1,numCompsTile);
compWeights = exp(-0.3*compCents);

% Make prior from components
priorTile = zeros(1,dx);

for ii = 1:numCompsTile
    
    priorTile = priorTile + compWeights(ii)*normpdf(supp,compCents(ii),compWidths(ii));
    
end


%% Plot both approaches

f1 = figure;
f1.Position = [100 100 1450 650];
hold on;

% Zero-mean
subplot(1,2,1);
hold on;

for ii = 1:numCompsZM
    if ii == numCompsZM
        p1 = plot(supp,wComp(ii)*normpdf(supp,muComp(ii),sigComp(ii)),'k','linewidth',2);
    else
        plot(supp,wComp(ii)*normpdf(supp,muComp(ii),sigComp(ii)),'k','linewidth',2);
    end
end
p2 = plot(supp,priorZM,'color',colorMat(2,:),'linewidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[suppLB suppUB],'xtick',[],'ytick',[]);
legend([p1,p2],{'Components','Prior'},'location','northwest','box','off');

% Fixed-components
subplot(1,2,2);
hold on;


for ii = 1:numCompsTile
    if ii == numCompsTile
        p3 = plot(supp,compWeights(ii)*normpdf(supp,compCents(ii),compWidths(ii)),'k','linewidth',2);
    else
        plot(supp,compWeights(ii)*normpdf(supp,compCents(ii),compWidths(ii)),'k','linewidth',2);
    end
end
p4 = plot(supp,priorTile,'color',colorMat(2,:),'linewidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[suppLB suppUB],'xtick',[],'ytick',[],'ylim',[0 4]);
legend([p3,p4],{'Components','Prior'},'location','northwest','box','off');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig9_MoGConstrainedFitting'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig9/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'MoGConstrainedFitting.svg']);
    
end

