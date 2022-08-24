%% Interactive script to show number of trials needed to accurately fit prior

clear all
close all

% Define number of datasets to simulate in interactive loop
numDatasets = 10;
countMin    = 5;
countMax    = 1000;
trialCounts = linspace(log(countMin),log(countMax),numDatasets);
trialCounts = round(exp(trialCounts));

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
mgrid     = -8:dx:8;                        % (just keep it fixed for paper figure)

% Define a prior
priortype = 3;  % valid options: 1 = gaussian, 2 = exp, 3 = cauchy, 4 = bimodal

switch priortype
    case 1, prior = normpdf(xgrid,0,2);   % gaussian
    case 2, prior = exp(-abs(xgrid))/2;   % exponential
    case 3, prior = (1./(1+xgrid.^2))/pi; % cauchy
    case 4, prior = normpdf(xgrid,-2,1)...
                    + normpdf(xgrid,2,1); % bimodal
end

prior = prior/sum(prior*dx); % normalize the prior to sum to 1
%------------------%

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

JSD = nan(numDatasets,1);

for ii = 1:numDatasets
%% Simulate data from the model 

% Randomly define a set of stimuli and randomly select a measurement value
% for each stimulus presentation
nsmps = trialCounts(ii);                                          % # of samples
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
[signsehat,priorhat,bwtshat,logliFinal,Mposthat] = fitEstimData_numerical(xdat,xhat,gbasis,xgrid,mgrid);

% inferred BLS estimate for each m value
BLSestimhat = Mposthat*xgrid*dx; 


%% Calculate divergence between true prior and estimated prior
JSD(ii) = getJSDiv(prior,priorhat);


%% Make plots

% True and inferred prior
%------------------% 
if ii == 1
    f1 = figure(1);
    f1.Position = [100 450 875 825];

    ymax  = 1.1*max([prior(:);priorhat(:)]);
    xlims = round(abs(xgrid(1)));
    ylims = round(abs(mgrid(1)));

    subplot(221) % linear scale
    hold on;
    p1 = plot(xgrid,prior,'r','linewidth',3);
    p2 = plot(xgrid,priorhat,'k--','linewidth',3);
    title('Prior'); box off;
    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],'xtick',linspace(-ylims,ylims,5),'ylim',[0 ymax]);
    legend([p1 p2],{'True prior', 'Inferred prior'},'AutoUpdate','off');
    xlabel('x');
    ylabel('p(x)');

    subplot(223) % log scale
    hold on;
    p3 = semilogy(xgrid,prior,'r','linewidth',3);
    p4 = semilogy(xgrid,priorhat,'k--','linewidth',3);
    title('Log-scale prior'); box off;
    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'yscale','log',...
        'xlim',ylims*[-1 1],'xtick',linspace(-ylims,ylims,5),'ylim',[1e-3 1]);
    legend([p3 p4],{'True prior', 'Inferred prior'},'AutoUpdate','off');
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
    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],...
        'xtick',linspace(-ylims,ylims,5),'ylim',ylims*[-1 1],'ytick',linspace(-ylims,ylims,5));
    legend('True BLS estimate', 'location', 'northwest','AutoUpdate','off');
    ylabel('Measurement (m)');

    subplot(224); % inferred posterior
    imagesc(xgrid,mgrid,Mposthat); axis xy;
    title('Inferred posterior');
    hold on;
    p5 = plot(BLSestim,mgrid,'k','linewidth',3);
    p6 = plot(BLSestimhat,mgrid,'c--','linewidth',3);
    % plot(xgrid,xgrid,'k--','linewidth',2);
    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',ylims*[-1 1],...
        'xtick',linspace(-ylims,ylims,5),'ylim',ylims*[-1 1],'ytick',linspace(-ylims,ylims,5));
    legend([p5 p6],{'True BLS', 'Inferred BLS'}, 'location', 'northwest','AutoUpdate','off');
    xlabel('Stimulus (x)');  ylabel('Measurement (m)');
else
    figure(1)
    subplot(221)
    plot(xgrid,priorhat,'k--','linewidth',3);
    subplot(223)
    semilogy(xgrid,priorhat,'k--','linewidth',3);
    subplot(224)
    plot(BLSestimhat,mgrid,'c--','linewidth',3);
end

% Fit quality (measured with Jensen-Shannon divergence)
%------------------% 
if ii == 1
    f2 = figure(2);
    f2.Position = [1000 600 650 600];
    hold on

    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 1000],'ylim',[1e-5 1],...
            'xtick',[0 10 100 500 1000],'yscale','log','xscale','log');
    xlabel('Number of trials');
    ylabel('Jensen-Shannon Divergence');
end

figure(2)
scatter(trialCounts(ii),JSD(ii),150,'k','filled');

if ii~=numDatasets
    disp(['Press a key to add ',num2str(trialCounts(ii+1)-trialCounts(ii)),' trials to dataset and refit!'])
    pause;
end

end