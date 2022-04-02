%% Figure 7 - Demo showing how to constrain number of parameters

clear all
close all

% Define numerical support
suppLB = -5;
suppUB = 5;
dx     = 300;
supp   = linspace(suppLB,suppUB,dx);

% Define color palette used for plotting
colorMat = colororder;

% Toggle on/off saving figures
saveOn = 1;


%% Restrict components to zero mean

% Define number of components to use
numCompsZM = 5;

% Define component parameters
wComp   = [0.24 0.08 0.16 0.23 0.29];
wComp   = wComp/sum(wComp);
muComp  = zeros(1,numCompsZM);
sigComp = [2.20 1.8 1.2 0.6 0.2];

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
f1.Position = [100 100 1250 650];
hold on;

% Zero-mean
subplot(1,2,1);
hold on;

for ii = 1:numCompsZM
    if ii == numCompsZM
        p1 = plot(supp,wComp(ii)*normpdf(supp,muComp(ii),sigComp(ii)),'k','linewidth',4);
    else
        plot(supp,wComp(ii)*normpdf(supp,muComp(ii),sigComp(ii)),'k','linewidth',4);
    end
end
p2 = plot(supp,priorZM,'color',colorMat(2,:),'linewidth',4);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',30,'xlim',[suppLB suppUB],'xtick',[],'ytick',[]);
legend([p1,p2],{'Components','Prior'},'location','northwest','box','off');
ylabel('Probability');
xlabel('Stimulus (x)');

% Fixed-components
subplot(1,2,2);
hold on;


for ii = 1:numCompsTile
    if ii == numCompsTile
        p3 = plot(supp,compWeights(ii)*normpdf(supp,compCents(ii),compWidths(ii)),'k','linewidth',4);
    else
        plot(supp,compWeights(ii)*normpdf(supp,compCents(ii),compWidths(ii)),'k','linewidth',4);
    end
end
p4 = plot(supp,priorTile,'color',colorMat(2,:),'linewidth',4);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',30,'xlim',[suppLB suppUB],'xtick',[],'ytick',[],'ylim',[0 4]);
legend([p3,p4],{'Components','Prior'},'location','northwest','box','off');
ylabel('Probability');
xlabel('Stimulus (x)');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig7_MoGConstrainedFitting'),filesep,'split');
    topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig7/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'Fig7_constraintExamples.svg']);
    
end

