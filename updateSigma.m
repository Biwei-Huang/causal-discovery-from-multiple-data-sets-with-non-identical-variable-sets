function sigma = updateSigma(mu, posqi, posCovEqi)
% update sigma

% the posterior of qi
posCovEqi = sum(posCovEqi,3);
posqi = sum(posqi,3);
sigma = posCovEqi./(posqi+eps^20)-mu.^2;
sigma = sigma.^0.5;