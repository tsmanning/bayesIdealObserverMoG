function [muProd,sigProd] = findGaussProd(mu1,sig1,mu2,sig2)

% Calculate the mean and sigma of the gaussian that results from
% multiplying two gaussians

muProd = (mu1.*sig2.^2 + mu2.*sig1.^2)./(sig2.^2 + sig1.^2);

sigProd = sqrt( ( (sig1.*sig2).^2 )./(sig2.^2 + sig1.^2) );

end

