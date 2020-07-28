function mu = updateMu(posqi, posMeanEqi)
% update mu

% the posterior of qi
posMeanEqi = sum(posMeanEqi,3);
posqi = sum(posqi,3);
mu = posMeanEqi./(posqi+eps^20);