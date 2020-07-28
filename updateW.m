function w = updateW(posqi)
% update w

% the posterior of qi
w = sum(posqi,3)/size(posqi,3);