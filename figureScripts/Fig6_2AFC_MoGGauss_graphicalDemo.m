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
saveOn = 1;


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
        % Don't recalculate values aready calculated
        pf1_plot(ii) = pf1;
        pf2_plot(ii) = pf2;
        pf3_plot(ii) = pf3;
    else
        % Define numerical supports
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
    topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
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

%% Helper function for plotting/calculating psychometric functions
function [p2AFC,db,f1] = get2AFCresponseP(x1,x2,sigL,wP,muP,sigP,plotOn)

%% Make sure parameters are column vectors
if ~iscolumn(wP)
    wP = wP';
end

if ~iscolumn(muP)
    muP = muP';
end

if ~iscolumn(sigP)
    sigP = sigP';
end


%% Define numerical support

% How many elements in 1D support?
numElem = 100;

% How many standard deviations to draw support out to?
facSD   = 4;

% Center joint likelihood over true stim and extend to 3SD
m1 = x1 + facSD*sigL(1)*linspace(-1,1,numElem);
m2 = x2 + facSD*sigL(2)*linspace(-1,1,numElem);

dx = [diff(m1(1:2)) diff(m2(1:2))];


%% Define helper functions
alphaFxn = @(gam,sig)           gam.^2 ./ (gam.^2 + sig.^2);
muFxn    = @(mu,gam,sig)        mu .* sig.^2 ./ (gam.^2 + sig.^2);
vFxn     = @(w,gam,mus,sig,muL) w .* 1./sqrt(gam.^2 + sig.^2) .* normcdf((muL - mus) ./ sqrt(gam.^2 + sig.^2));
wFxn     = @(w,gam,mus,sig,muL) vFxn(w,gam,mus,sig,muL)./sum(vFxn(w,gam,mus,sig,muL),2);

postFxn  = @(supp,w,gam,nus,sig,m) ...
            wFxn(w,gam,nus,sig,m)'*normpdf(supp,alphaFxn(gam,sig)*m + muFxn(nus,gam,sig),sqrt(alphaFxn(gam,sig))*sig);


%% Get p(x2>x1|m1,m2) for each measurement pair in grid

% Stim 1
for ii = 1:numElem

    % Stim 2
    for jj = 1:numElem

        % Make a grid of STIMULUS values about current measurements
        thisM1 = m1(ii);
        thisM2 = m2(jj);
        
        maxSig1 = max(sqrt(alphaFxn(sigP,sigL(1)))*sigL(1));
        maxSig2 = max(sqrt(alphaFxn(sigP,sigL(2)))*sigL(2));
        
        theseX1 = thisM1 + facSD*maxSig1*linspace(-1,1,numElem);
        theseX2 = thisM2 + facSD*maxSig2*linspace(-1,1,numElem);
        [xM1,xM2] = meshgrid(theseX1,theseX2);
        
        % Compute posterior over grid for this set of measurements
        thisPost1 = postFxn(theseX1,wP,sigP,muP,sigL(1),thisM1);
        thisPost2 = postFxn(theseX2,wP,sigP,muP,sigL(2),thisM2);
        
        this2Dpost = thisPost2'*thisPost1;
        this2Dpost = this2Dpost/sum(this2Dpost(:));
        
        % Sum probability over unity line to get p(x2>x1|m1,m2)
        thisMask = xM2>xM1;
        maskedPost = this2Dpost.*thisMask;
        
        post2AFC(jj,ii) = sum(maskedPost(:));
        
    end
    
end


%% Define decision boundary and mask [i.e. p(x2>x1|m1,m2) > 0.5]

% Want to also include half the probability from bins falling on the
% boundary

% decMask units: columns-m2 min to m2 max, rows- m1 min to m2 max

decMask = post2AFC > 0.5;

db1     = diff(decMask,[],1);
db2     = diff(decMask,[],2);
dbMask  = or([zeros(1,numElem);db1]<0,[zeros(numElem,1) db2]<0) == 1;

[dbyInd,dbxInd] = find(dbMask);
db = [m1(dbxInd);m2(dbyInd)];


%% Get p(x2>x1|x1,x2) by summing joint likelihood above the decision boundary

% Make joint likelihood
like2D = normpdf(m2,x2,sigL(2))'*normpdf(m1,x1,sigL(1));
like2D = like2D/sum(like2D(:));

maskedP1 = like2D.*decMask;
maskedP2 = like2D.*dbMask;
p2AFC = sum(maskedP1(:)) + 0.5*sum(maskedP2(:));


%% Plot if desired

if plotOn
   % Define gamma'd colormap to use for plots
   cmap = repmat(linspace(0,1,256)'.^1.0,[1 3]);
    
   lb = max([m1(1) m2(1)]);
   ub = min([m1(end) m2(end)]);
    
   f1 = figure;
   f1.Position = [100 100 650 650];
   hold on;
   
   imagesc(m1,m2,like2D); axis xy; axis image;
   plot(db(1,:),db(2,:),'w','linewidth',5);
   set(gca,'plotboxaspectratio',[1 1 1],'xtick',[],'ytick',[],'xlim',[lb ub],'ylim',[lb ub],'fontsize',20);
%    set(gca,'xlim',[m1(1) m1(end)],'ylim',[m2(1) m2(end)],'fontsize',25);
   xlabel('Stimulus 1 (x_{1})');
   ylabel('Stimulus 2 (x_{2})');
   text(db(1,1),db(2,28),'Yes','color',[1 1 1],'fontsize',40);
   text(db(1,16),db(2,8),'No','color',[1 1 1],'fontsize',40);
   colormap(cmap)

else
    
    f1 = [];
   
end

end