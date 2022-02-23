%% Figure 3: Simulate an observer's estimate distribution from many trials 

clear all
close all

% Define numerical support
suppLB  = -7;                                   % lower bound
suppUB  = 7;                                    % upper bound
numBins = 300;                                  % number of bins
supp    = linspace(suppLB,suppUB,numBins);      % stimulus support
dx      = diff(supp(1:2));                      % precision

% Grab matlab color palette
colorMat = colororder;                          

% Toggle on/off saving figures
saveOn = 0;


%% Define zero-mean MoG prior

% Define number of components to use
numCompsZM = 5;

% Define component parameters
pW   = [0.2451 0.1680 0.0278 0.2847 0.2745];
pW   = pW/sum(pW);
pNu  = zeros(1,numCompsZM);
pGam = [0.6407 2.1001 2.7232 1.3693 2.1668];

% Make prior from components
prior = getMoGPrior(pGam,pNu,pW,supp);


%% Generate a set of likelihoods

% Define number of likelihoods
numLike = 5;

% Define expected value of likelihood and standard deviation
likeMean  = 3;
likeSig   = 1;
 
% Select some evenly spaced means for individual trial likelihoods
likeNoise = norminv(linspace(0.25,0.75,numLike),likeMean,likeSig);

likeMat = nan(numLike,numBins);

for ii = 1:numLike
    
    likeMat(ii,:) = normpdf(supp,likeNoise(ii),likeSig); 
    
end

% Likelihood centered on EV
expectLike = normpdf(supp,likeMean,likeSig);


%% Generate set of posteriors

% Get normalized product of prior and likelihoods

for ii = 1:numLike
    
    postMat(ii,:) = prior.*likeMat(ii,:); 
    postMat(ii,:) = postMat(ii,:)/(sum(postMat(ii,:))*dx);
    
end

% Posterior for expected value of stimulus
expectPost = prior.*expectLike;
expectPost = expectPost/(sum(expectPost)*dx);


%% Generate estimate distribution

% Generate a set of measurements
nMeas = 100000;
meas  = randn(nMeas,1) + likeMean;

% Get Bayes estimates (center of mass)
for ii = 1:nMeas
    
    thisPost = prior.*normpdf(supp,meas(ii),likeSig);
    thisPost = thisPost/sum(thisPost);
    xHat(ii) = thisPost*supp';
    
end

% Make a histogram of the estimates and normalize to PMF
[estDist,edges] = histcounts(xHat,35);
binCents        = edges(1:end-1) + diff(edges(1:2));
estDist         = estDist/(sum(estDist)*diff(edges(1:2)));


%% Plot

maxPost = 1.05*max([postMat(:)']);

% Prior + like
f1 = figure;
f1.Position = [100 100 650 650];
hold on

plot([supp(1) supp(end)],[0 0],'k','linewidth',2);
p1 = plot(supp,prior,'color',[0.9 0 0],'linewidth',6);

p2 = plot(supp,likeMat(1,:),'color',[0.7 0.7 0.7],'linewidth',6);

for ii = 2:numLike
   
    plot(supp,likeMat(ii,:),'color',[0.7 0.7 0.7],'linewidth',6);
    
end

p3 = plot(supp,expectLike,'color',[0 0 0],'linewidth',6);
p3a = arrow([likeMean -0.125*maxPost],[likeMean -0.025*maxPost],...
        'color',[0 0 0],'linewidth',5);
    
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[supp(1) supp(end)],'ylim',[-0.125*maxPost maxPost],'ytick',[]);
xlabel('Stimulus');
ylabel('Probability');
legend([p1,p2,p3,p3a],{'Prior','Likelihood','Likelihood (EV)','True stimulus'},'location','northwest');

% Prior + posts
f2 = figure;
f2.Position = [800 100 735 650];
hold on

plot([supp(1) supp(end)],[0 0],'k','linewidth',2);
p4 = plot(supp,prior,'color',[0.9 0 0],'linewidth',6);

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[supp(1) supp(end)],'ylim',[-0.125*maxPost maxPost],'ytick',[]);

p5 = plot(supp,postMat(1,:),'color',[0.6 0.6 1],'linewidth',6);
p5a = arrow([supp*postMat(1,:)'/sum(postMat(1,:)) 0.125*maxPost],[supp*postMat(1,:)'/sum(postMat(1,:)) 0.025*maxPost],...
        'color',[0.6 0.6 1],'linewidth',5);

for ii = 2:numLike
   
    plot(supp,postMat(ii,:),'color',[0.6 0.6 1],'linewidth',6);
    arrow([supp*postMat(ii,:)'/sum(postMat(ii,:)) 0.125*maxPost],[supp*postMat(ii,:)'/sum(postMat(ii,:)) 0.025*maxPost],...
        'color',[0.6 0.6 1],'linewidth',5);

end

p6 = plot(supp,expectPost,'color',[0 0 0.85],'linewidth',6);
arrow([supp*expectPost'/sum(expectPost) 0.125*maxPost],[supp*expectPost'/sum(expectPost) 0.025*maxPost],...
    'color',[0 0 0.85],'linewidth',5);

xlabel('Stimulus'); ylabel('Probability');
legend([p4,p5,p6],{'Prior','Posterior','Posterior (EV)'},'location','northwest');


% Estimate distribution
maxEstD = maxPost;

f3 = figure;
f3.Position = [1000 100 725 675];
hold on;

plot([supp(1) supp(end)],[0 0],'k','linewidth',2);
p8 = bar(binCents,estDist,'EdgeColor',[0 0 0],'FaceColor',[0.6 0.6 1],'barwidth',1);

p9 = arrow([supp*postMat(1,:)'/sum(postMat(1,:)) -0.125*maxEstD],[supp*postMat(1,:)'/sum(postMat(1,:)) -0.025*maxEstD],...
        'color',[0.6 0.6 1],'linewidth',5);
for ii = 2:numLike
   
    arrow([supp*postMat(ii,:)'/sum(postMat(ii,:)) -0.125*maxEstD],[supp*postMat(ii,:)'/sum(postMat(ii,:)) -0.025*maxEstD],...
        'color',[0.6 0.6 1],'linewidth',5);

end
p10 = arrow([supp*expectPost'/sum(expectPost) -0.125*maxEstD],[supp*expectPost'/sum(expectPost) -0.025*maxEstD],...
    'color',[0 0 0.85],'linewidth',5);

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[supp(1) supp(end)],'ytick',[],'ylim',[-0.125*maxEstD maxEstD]);
xlabel('Stimulus');
ylabel('p(xHat|x)');
legend([p8,p9,p10],{'Estimate Distribution','Estimate','Estimate (EV)'},'location','northwest');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig3_EstimateDistribution'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig3/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'priorLike.svg']);
    saveas(f2,[sDir,'priorPost.svg']);
    saveas(f3,[sDir,'priorEstDist.svg']);
    
end
