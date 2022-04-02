%% Fig 8 - Fit an observer model with a constrained mixture of Gaussians prior with estimate data
% CAUCHY PRIOR

close all
clear all

% Define observer/stimulus parameters
%------------------% 
% Define observation noise 
lSig      = 1.2;                            % standard deviation 

% Define numerical support for measurements and stimulus values
xtestrnge = 4*[-1 1];                       % range of test stimuli to consider
dx        = .1;                             % bin size
mrnge     = xtestrnge + 3*lSig*[-1 1];      % range of measurements to consider
xrnge     = mrnge + 2*lSig*[-1 1];          % range of posterior over x
xgrid     = (xrnge(1)+dx/2:dx:xrnge(2))';   % stimulus grid
% mgrid = (mrnge(1)+dx/2:dx:mrnge(2))';     % internal measurement grid
mgrid     = -8:dx:8;                        % (just keep it fixed for paper figure)

% Define a prior
priortype = 3;  % valid options: 1 = gaussian, 2 = exp, 3 = cauchy

switch priortype
    case 1, prior = normpdf(xgrid,0,2);   % gaussian
    case 2, prior = exp(-abs(xgrid))/2;   % exponential
    case 3, prior = (1./(1+xgrid.^2))/pi; % cauchy
end

prior = prior/sum(prior*dx); % normalize the prior to sum to 1
%------------------%

% Toggle on/off saving figures
saveOn = 0;


%% Compute posterior and BLS estimate for all possible measurements m

% Make a grid of stim and measurement values
[xx,mm]  = meshgrid(xgrid,mgrid); 

% Calculate likelihoods for each stimulus value
ppli     = normpdf(mm,xx,lSig);

% Calulate posteriors for each likelihood
ppost    = ppli.*prior';             % unnormalized posterior
ppost    = ppost./(sum(ppost,2)*dx); % normalized posterior

% Calculate BLS estimate for each posterior
BLSestim = ppost*xgrid*dx;           % BLS estimate for each m value


%% Simulate data from the model 

% Randomly define a set of stimuli and randomly select a measurement value
% for each stimulus presentation
nsmps = 1000;                                          % # of samples
xdat  = rand(nsmps,1)*diff(xtestrnge) + xtestrnge(1);  % stimuli
mdat  = xdat + randn(nsmps,1)*lSig;                    % observer's noisy measurements

% Generate a set of observer estimates (assuming the observer makes a BLS
% estimate) by interpolating between values in the m-xHat lookup table
xhat = interp1(mgrid,BLSestim,mdat,'linear','extrap');


%% Set up non-fitted parameters of constrained MoG prior

% Define Gaussian basis set 
nbasis  = 6;                                 % number of Gaussian basis functions for prior 
bctrs   = zeros(1,nbasis);                   % prior means
sigshft = 2;                                 % shift the location of first sigma
bsigs   = (2.^(0-sigshft:nbasis-1-sigshft)); % prior standard deviations

% Make basis set
basisFun = @(x,mus,sigs)(exp(-0.5*(x-mus).^2./sigs.^2)./(sqrt(2*pi)*sigs));

gbasis   = basisFun(xgrid,bctrs,bsigs); % basis centers
gbasis   = gbasis./(dx*sum(gbasis));    % normalize so each sums to 1


%% Fit the model to the synthetic data set

% Use maximum likelihood estimation to determine best fitting weights for
% basis functions of MoG prior given the dataset
[signsehat,priorhat,bwtshat,logliFinal,Mposthat] = fitBLSobserverModel_estimdata(xdat,xhat,gbasis,xgrid,mgrid);

% inferred BLS estimate for each m value
BLSestimhat = Mposthat*xgrid*dx; 


%% Make plots

% True and inferred prior
%------------------% 
f1 = figure;
f1.Position = [100 100 875 825];

ymax  = 1.1*max([prior(:);priorhat(:)]);
xlims = round(abs(xgrid(1)));
ylims = round(abs(mgrid(1)));

subplot(221) % linear scale
hold on;
plot(xgrid,prior,'r','linewidth',3);
plot(xgrid,priorhat,'k--','linewidth',3);
title('Prior'); box off;
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],'xtick',linspace(-ylims,ylims,5),'ylim',[0 ymax]);
legend('True prior', 'Inferred prior');
xlabel('x'); 
ylabel('p(x)');

subplot(223) % log scale
hold on;
semilogy(xgrid,prior,'r','linewidth',3);
semilogy(xgrid,priorhat,'k--','linewidth',3);
title('Log-scale prior'); box off;
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'yscale','log',...
    'xlim',ylims*[-1 1],'xtick',linspace(-ylims,ylims,5),'ylim',[1e-3 1]);
legend('True prior', 'Inferred prior');
xlabel('x'); 
ylabel('p(x)');

% Plot true and inferred posterior & BLS estimates
%------------------% 
subplot(222);  % true posterior
imagesc(xgrid,mgrid,ppost); axis xy;
title('True posterior');
hold on;
plot(BLSestim,mgrid,'k','linewidth',3); 
% plot(xgrid,xgrid,'k--','linewidth',2);
hold off;
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],...
    'xtick',linspace(-ylims,ylims,5),'ylim',ylims*[-1 1],'ytick',linspace(-ylims,ylims,5));
legend('True BLS estimate', 'location', 'northwest');
ylabel('Measurement (m)');

subplot(224); % inferred posterior
imagesc(xgrid,mgrid,Mposthat); axis xy;
title('Inferred posterior');
hold on;
plot(BLSestim,mgrid,'k',BLSestimhat,mgrid,'c--','linewidth',3); 
% plot(xgrid,xgrid,'k--','linewidth',2);
hold off;
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],...
    'xtick',linspace(-ylims,ylims,5),'ylim',ylims*[-1 1],'ytick',linspace(-ylims,ylims,5));
legend('True BLS', 'Inferred BLS', 'location', 'northwest');
xlabel('Stimulus (x)');  ylabel('Measurement (m)');


%% Save figures

if saveOn
    
    splPath = regexp(which('Fig8_MoGtoNonGauss'),filesep,'split');
    topDir  = [fullfile(splPath{1:numel(splPath)-1}),filesep];
    sDir = [topDir,'figuresImgs/fig8/'];
    
    if ~isfolder(sDir)
        mkdir(sDir)
    end
    
    saveas(f1,[sDir,'MoGtoNonGauss.svg']);
    
end

