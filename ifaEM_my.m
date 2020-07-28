function [A, mu, sigma, w, loglAll] = ifaEM_my(X, IDX0, numGauss, parsEM, filename)

% EM algorithm to learn the mixing matrix over the intergrated set of variables
% Inputs:
%       X: data points
%       IDX0: each column denotes the indices of variables in each data set; 
%       numGauss: number of Gaussian for each component
%       parsEM: parameters for EM
%       filename: the file name to save the results

% Outputs:
%       A: mixing matrix
%       mu: mixture of gaussian means
%       sigma: mixture of gaussian stds for all components
%       w: mixture of gaussian weights
%       loglAll: log-likelihood


M = length(X); % number of data sets
n=size(X{1},1); % total number of integrated variables
thres = parsEM.thres;

% initialize parameters
A = parsEM.A;
for m = 1:M
    mu{m} = parsEM.mu;
    sigma{m} = parsEM.sigma;
    w{m} = parsEM.w;
    Dim{m} = size(w{m},1);
end
Lambda = parsEM.noise*eye(n);

iter = 0;
loglAll = [];
logl_prev = inf;
logl2 = inf;
while iter == 0 || iter < parsEM.maxIter && abs(logl_prev-logl2) >= abs(logl_prev)*thres
    % M step
    % update A
    logl_prev = logl2;
    if iter > 0
%         A = (X*posMeanE')/(sum(posCovE,3));
        
        %%%%%%%%%%%%%%%%
        tmp1 = zeros(n,n);
        for m = 1:M % total number of data sets
           tmp1 = tmp1 + X{m}*posMeanE{m}';
        end
        A = zeros(n,n);
        for i = 1:n % total number of variables
           tmp2 = zeros(n,n);
           idx = find(IDX0(i,:)==1); % the index of the datasets which have the i-th variable
           for j = idx
              tmp2 = tmp2 + sum(posCovE{j},3);
           end
           A(i,:) = tmp1(i,:)/(tmp2+ 1e-4*eye(n));
        end
        
        
        if parsEM.updatePrior == 1
            for m = 1:M  % total number of data sets
                posqi{m} = marginalPosQ(posq{m}, Dim{m}, 0, numGauss, qidMat{m});
                % update mu
                if mu{m}~=zeros(Dim{m},numGauss)
                    mu{m} = updateMu(posqi{m}, posMeanEqi{m});
                end
                % update sigma
                sigma{m} = updateSigma(mu{m}, posqi{m}, posCovEqi{m});
                
                % update w
                w{m} = updateW(posqi{m});
            end
        end
    end
    
    % E step
    for m = 1:M % total number of data sets
        qidMat{m} = qidMatrix(numGauss, Dim{m}, 0);
        priorq{m} = priorLatentQ(w{m}, Dim{m}, 0, numGauss, qidMat{m});
        
        % posterior of q and e
        I{m} = diag(IDX0(:,m));
        [posq{m}, marginalx{m}] = evaluateQ(X{m}, I{m}*A, mu{m}, sigma{m}, Dim{m}, 0, numGauss, qidMat{m}, priorq{m},Lambda);
        [posMeanE{m}, posCovE{m}, posMeanEqi{m}, posCovEqi{m}] = evaluateE_my(X{m}, posq{m}, I{m}*A, Lambda, ...
            mu{m}, sigma{m}, Dim{m}, 0, numGauss, qidMat{m});
        logl(m) = -sum(log(marginalx{m} + 1e-300*ones(length(marginalx{m}),1)));
    end
    logl2 = sum(logl);
    fprintf('iter%d: negative loglik: %.10f\n', iter, logl2);
%     A
    iter = iter + 1;
    loglAll = [loglAll, logl2];
    
    save(filename,'IDX0','A');
end